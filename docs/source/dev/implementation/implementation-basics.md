# The Basics

Vistle is designed to be easily extensible through its modular structure. This tutorial will show you how to write your own Vistle module and how to incorporate it into the software, so you can use it right away.

## How to write a general module
- skeleton: compute function..., input and output ports, how does it translate to GUI
- this is the bare minimum, also talk about other classes in vistle::Module!

`MyModule.h`
```cpp
#ifndef MYMODULE_H
#define MYMODULE_H

#include <vistle/module/module.h>

// this is a test module for the developer guide
class MyModule: public vistle::Module {
public:
    MyModule(const std::string &name, int moduleID, mpi::communicator comm);
    ~MyModule();

private:
    bool compute(const std::shared_ptr<vistle::BlockTask> &task) const override;
};

#endif // MYMODULE_H
```

`MyModule.cpp`:

```cpp
#include <vistle/alg/objalg.h>
#include <vistle/core/uniformgrid.h>

#include "MyModule.h"

MODULE_MAIN(MyModule)

using namespace vistle;

MyModule::MyModule(const std::string &name, int moduleID, mpi::communicator comm): Module(name, moduleID, comm)
{
    createInputPort("data_in", "input grid with mapped data");
    createOutputPort("grid_out", "output grid");
}

MyModule::~MyModule()
{}

bool MyModule::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
    auto input = task->expect<Vec<Scalar, 3>>("data_in");
    if (!input) {
        sendError("This module only supports three-dimensional vector data fields!");
        return true;
    }

    auto container = splitContainerObject(input);
    auto grid = container.geometry;
    if (!grid) {
        sendError("Found no grid at the input port!");
        return true;
    }

    UniformGrid::ptr outputGrid(new UniformGrid(1, 1, 1));
    Vec<Scalar, 1>::ptr outputData(new Vec<Scalar, 1>(1));

    outputGrid->copyAttributes(grid);
    updateMeta(outputGrid);

    outputData->setGrid(outputGrid);
    updateMeta(outputData);

    task->addObject("grid_out", outputGrid);

    return true;
}
```

`CMakeLists.txt`:
```cmake
add_module(MyModule "Test module" MyModule.h MyModule.cpp)
```

## How to write a reader module
- create a MyReader test module similar to MyModule above...
- this is the bare minimum also talk about other classes in vistle::Reader!

`MyReader.h`:

```cpp
#ifndef MYREADER_H
#define MYREADER_H

#include <vistle/module/reader.h>

class MyReader: public vistle::Reader {
public:
    MyReader(const std::string &name, int moduleID, mpi::communicator comm);
    ~MyReader() override;

    bool read(vistle::Reader::Token &token, int timestep = -1, int block = -1) override;
};

#endif //MYREADER_H
```

`MyReader.cpp`:

```cpp
#include "MyReader.h"

MODULE_MAIN(MyReader)

using vistle::Reader;

MyReader::MyReader(const std::string &name, int moduleID, mpi::communicator comm): Reader(name, moduleID, comm)
{
    createOutputPort("grid_out", "output grid");
}

MyReader::~MyReader()
{}

bool MyReader::read(vistle::Reader::Token &token, int timestep, int block)
{
    return true;
}
```

`CMakeLists.txt`
```
add_module(MyReader "test read module" MyReader.h MyReader.cpp)
```

## How to add your module to Vistle

- create and add everything to `module/test/MyModule` folder (--> explain different module categories)
- `MyReader` goes to `module/read/MyReader`

`module/test/CMakeLists.txt`:

```cmake
# modules mainly suitable for testing and developing
set(VISTLE_MODULE_CATEGORY "Test")

add_subdirectory(ClipVtkm)
# ...
add_subdirectory(MpiInfo)
add_subdirectory(MyModule) # add this line to include our test module in vistle
add_subdirectory(ObjectStatistics)
# ...

```
## How to add your module to the Vistle release
- maybe explain how to do a pull request to add module to vistle?
- how to create documentation for the module (link to page in documentation that does not exist yet)
