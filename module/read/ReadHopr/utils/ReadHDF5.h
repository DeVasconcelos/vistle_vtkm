#ifndef READHOPR_READHDF5_H
#define READHOPR_READHDF5_H

#include <algorithm>
#include <cassert>
#include <hdf5.h>
#include <iostream>
#include <string>
#include <vector>

#include <vistle/core/parameter.h>

template<typename T>
struct H5Dataset {
    std::vector<T> vector;
    std::vector<hsize_t> dimension;
};

//TODO: strings must be handled differently (as they can be variable or fixed size, see 'readH5Attribute<std::string>')
template<typename T>
H5Dataset<T> readH5Dataset(hid_t fileId, const char *datasetName)
{
    H5Dataset<T> result;
    auto datasetId = H5Dopen(fileId, datasetName, H5P_DEFAULT);
    if (datasetId < 0) {
        std::cerr << "Could not open dataset " << datasetName << "!" << std::endl;
        H5Dclose(datasetId);
        return result;
    }
    auto spaceId = H5Dget_space(datasetId);
    std::vector<hsize_t> shape(H5Sget_simple_extent_ndims(spaceId));
    H5Sget_simple_extent_dims(spaceId, shape.data(), nullptr);

    for (hsize_t dim: shape) {
        result.dimension.push_back(dim);
    }

    auto totalSize = std::accumulate(result.dimension.begin(), result.dimension.end(), 1, std::multiplies<hsize_t>());

    auto dtypeId = H5Dget_type(datasetId);
    result.vector.resize(totalSize);
    if (H5Dread(datasetId, dtypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, result.vector.data()) < 0) {
        std::cerr << "An error occurred when trying to read in " << datasetName << " of type " << dtypeId << "."
                  << std::endl;
        result.vector.clear();
    }

    //TODO: find a way to not always have to close all ids manually
    H5Tclose(dtypeId);
    H5Sclose(spaceId);
    H5Dclose(datasetId);
    return result;
}

//TODO: use struct like H5Dataset here, too!
template<typename T>
std::vector<T> readH5Attribute(hid_t fileId, const char *attrName)
{
    std::vector<T> result;
    hsize_t attrSize = 1;
    if (H5Aexists(fileId, attrName)) {
        auto attrId = H5Aopen(fileId, attrName, H5P_DEFAULT);
        if (attrId > -1) {
            auto spaceId = H5Aget_space(attrId);
            if (spaceId > 1) {
                std::vector<hsize_t> shape(H5Sget_simple_extent_ndims(spaceId));
                H5Sget_simple_extent_dims(spaceId, shape.data(), nullptr);

                for (hsize_t dim: shape) {
                    attrSize *= dim;
                }

                result.resize(attrSize);
                if (H5Aread(attrId, H5Aget_type(attrId), result.data()) < 0) {
                    std::cerr << "An error occurred when trying to read the attribute '" << attrName
                              << "' with H5Aread!" << std::endl;
                    H5Sclose(spaceId);
                    H5Aclose(attrId);
                    return std::vector<T>();
                }
            }
            H5Sclose(spaceId);
        } else {
            std::cerr << "An error occurred when trying to open the attribute '" << attrName << "'with H5Aopen."
                      << std::endl;
            return result;
        }
        H5Aclose(attrId);
    } else {
        std::cerr << "The attribute '" << std::to_string(attrName) << "' does not exist in the given file! (H5Aexists)"
                  << std::endl;
    }
    return result;
}

// HDF5 string attributes can be of fixed or variable length, so they have to be read in differently
// than attributes of other types.
template<>
std::vector<std::string> readH5Attribute<std::string>(hid_t fileId, const char *attrName)
{
    std::vector<std::string> result;
    hsize_t attrSize = 1;
    if (H5Aexists(fileId, attrName)) {
        auto attrId = H5Aopen(fileId, attrName, H5P_DEFAULT);
        if (attrId > -1) {
            auto spaceId = H5Aget_space(attrId);
            if (spaceId > 1) {
                std::vector<hsize_t> shape(H5Sget_simple_extent_ndims(spaceId));
                H5Sget_simple_extent_dims(spaceId, shape.data(), nullptr);

                for (hsize_t dim: shape) {
                    attrSize *= dim;
                }

                auto dtypeId = H5Aget_type(attrId);
                assert(H5Tget_class(dtypeId) == H5T_STRING);

                if (H5Tis_variable_str(dtypeId)) {
                    std::vector<char *> tmpResult(attrSize);

                    if (H5Aread(attrId, dtypeId, tmpResult.data()) < 0) {
                        std::cerr << "An error occurred when trying to read the attribute '" << attrName
                                  << "' with H5Aread!" << std::endl;
                        H5Sclose(spaceId);
                        H5Aclose(attrId);
                        return result;
                    }

                    for (hsize_t i = 0; i < attrSize; ++i) {
                        result.push_back(std::string(tmpResult[i]));
                        free(tmpResult[i]);
                    }
                } else {
                    size_t strSize = H5Tget_size(dtypeId);
                    std::vector<char> temp_data(attrSize * strSize);

                    if (H5Aread(attrId, dtypeId, temp_data.data()) < 0) {
                        std::cerr << "An error occurred when trying to read the attribute '" << attrName
                                  << "' with H5Aread!" << std::endl;
                        H5Sclose(spaceId);
                        H5Aclose(attrId);
                        return result;
                    }

                    for (hsize_t i = 0; i < attrSize; ++i) {
                        auto rawStr = std::string(&temp_data[i * strSize], strSize);
                        // remove whitespaces from string
                        rawStr.erase(std::remove_if(rawStr.begin(), rawStr.end(), ::isspace), rawStr.end());
                        result.push_back(rawStr);
                    }
                }
            }
            H5Sclose(spaceId);
        } else {
            std::cerr << "An error occurred when trying to open the attribute '" << attrName << "'with H5Aopen."
                      << std::endl;
            return result;
        }
        H5Aclose(attrId);
    } else {
        std::cerr << "The attribute '" << std::to_string(attrName) << "' does not exist in the given file! (H5Aexists)"
                  << std::endl;
    }
    return result;
}


#endif //READHOPR_READHDF5_H
