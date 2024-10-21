# The Basics

Vistle is designed to be easily extensible through its modular structure. This tutorial will show you how to write your own Vistle module and how to incorporate it into the software, so you can use it right away.

It will explain the following topics:
- [How to write a compute module](#how-to-write-a-compute-module)
- [How to write a reader module](#how-to-write-a-reader-module)
- [How to add your module to Vistle](#how-to-add-your-module-to-vistle)
- [How to add your module to the Vistle release](#how-to-add-your-module-to-the-vistle-release)

## How to write a compute module
<!-- TODO: add a section about the other functions in vistle::module (changeParameter, prepare, reduce, examine, ...) -->

The most common Vistle modules are the compute modules. They receive data from other modules through their input port(s), run computations on the data, and finally return the result(s) through their output port(s).

The following example shows the minimum code necessary to develop a Vistle module. It will create a module named **MyModule**. By convention, all source code files for one module are stored in a folder named after the module. 

### Overview
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

To create a new Vistle module, create a class which inherits from the `vistle::Module` class. The name of that class is the module's name and will also be displayed in the GUI later on.
The input and output ports as well as the parameter of the module are defined in the class's constructor. The computations done by the module are defined in the `compute` method. This method is called whenever the Vistle pipeline is executed or you double-click on the module.

### The main function

It's important to add the following line to the top of your module's source file. It automatically creates the main function of the module integrating it into Vistle correctly.
```cpp
MODULE_MAIN(MyModule)
```

### Defining ports
 The input and output ports are created by calling the `createInputPort` and `createOutputPort` functions in the module's constructor. Two strings, the name of the port (which can be used to reference it in the code later on) and a short description, must be passed to the function, as shown to the following example:
```cpp
MyModule::MyModule(const std::string &name, int moduleID, mpi::communicator comm): Module(name, moduleID, comm)
{
    // ports
    createInputPort("data_in", "input grid with mapped data");
    createOutputPort("grid_out", "output grid");
}
```

This code snippet adds an input and an output port to the module in the GUI:

<p align="center">
    <img src="myModulePorts.png" alt="MyModulePorts" width=120>
</p>
Hovering over either port will show the port's name and description.

### Working with ports
The following code snippet is an example of how to work with the data provided by the input port in the `compute` function.

```cpp
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

    return true;
}
```
The code checks if the data field attached to the input grid is a three-dimensional vector and if the input grid is valid. If not, an error message is printed using the `sendError` method. The error message will be shown in the Vistle console as well as in the parameter window. 

<p align="center">
    <img src="myModuleError.png" alt="MyModuleErrorMessage" >
</p>


The following code snippet shows how to add a uniform grid and a scalar data field to the module's output port:
```cpp
bool MyModule::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
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

### Defining module parameters
Module parameters can be defined by using the `addIntParameter`, `addFloatParameter`, `addStringParameter`, `addVectorParameter`, `addIntVectorParameter` or `addStringVectorParameter` functions in the module's constructor. All of these functions expect a string with the parameter's name, a string with the parameter's description and a default value. Some functions expect to be passed additional parameters or have overloaded function definitions.

In the following example code, we show how to add the most commonly used module parameters to **MyModule**. Note that class members for each of the module's parameters have been added to the header file (which is not shown here).

```cpp
DEFINE_ENUM_WITH_STRING_CONVERSIONS(Option, (Option1)(Option2)(Option3))

MyModule::MyModule(const std::string &name, int moduleID, mpi::communicator comm): Module(name, moduleID, comm)
{
    // parameters
    m_boolean = addIntParameter("boolean", "a boolean parameter", false, Parameter::Boolean);

    m_scalar = addFloatParameter("scalar", "a float parameter", 1.0);
    setParameterRange(m_scalar, 0.0, 10.0);

    m_vector = addVectorParameter("vector", "a vector parameter", ParamVector(1.0, 2.0, 3.0));

    m_choice = addIntParameter("choice", "choose one of the options", Option1, Parameter::Choice);
    V_ENUM_SET_CHOICES(m_choice, Option);
}
```
This code snippet created the following parameters which can be modified in the GUI's parameter window:
<p align="center">
    <img src="myModuleParameters.png" alt="MyModuleParameters" >
</p>

The first parameter is a boolean parameter with the default value false. Passing `Parameter::Boolean` to the `addIntParameters` function (as shown in the example), will create a checkbox that the user can toggle in the GUI.

The second parameter is a float value which the user can set in the parameter window. If desired, a minimum and maximum value can be defined for any arithmetic parameter using the `setParameterRange` method. The `scalar` parameter in this example only allows float values between 0 and 10.

The third parameter is a vector parameter which is created similarly to the scalar parameters. The user can define each component separately.

Finally, the fourth parameter is a choice parameter which creates a drop-down menu in the GUI which allows the user to set one of the available options. To create a choice parameter:
1. pass `Parameter::Choice` to the `addIntParameters` function
2. create an enum including all possible options using the `DEFINE_ENUM_WITH_STRING_CONVERSIONS` macro. Note that the names you define here will be displayed in the drop-down menu later on.
3. add the enum to the choice parameter with the `V_ENUM_SET_CHOICES` function.

### Working with module parameters

You can read the value of a module parameter using the `getValue()` function.

```cpp
bool MyModule::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
    // parameters
    bool boolean = m_boolean->getValue();
    Float scalar = m_scalar->getValue();
    ParamVector vector = m_vector->getValue();
    Option choice = (Option)m_choice->getValue();

    return true;
}
```

### Code
<details>
<summary> Click on the arrow on the left to view the complete MyModule header file. </summary>

`module/test/MyModule/MyModule.h`:

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

    vistle::IntParameter *m_boolean;
    vistle::FloatParameter *m_scalar;
    vistle::VectorParameter *m_vector;
    vistle::IntParameter *m_choice;
};

#endif // MYMODULE_H
```
</details>

<details>
<summary> Click on the arrow on the left to view the complete MyModule source file.</summary>

`module/test/MyModule/MyModule.cpp`

```cpp
#include <vistle/alg/objalg.h>
#include <vistle/core/uniformgrid.h>

#include "MyModule.h"

MODULE_MAIN(MyModule)

using namespace vistle;

DEFINE_ENUM_WITH_STRING_CONVERSIONS(Option, (Option1)(Option2)(Option3))

MyModule::MyModule(const std::string &name, int moduleID, mpi::communicator comm): Module(name, moduleID, comm)
{
    // ports
    createInputPort("data_in", "input grid with mapped data");
    createOutputPort("grid_out", "output grid");

    // parameters
    m_boolean = addIntParameter("boolean", "a boolean parameter", false, Parameter::Boolean);

    m_scalar = addFloatParameter("scalar", "a float parameter", 1.0);
    setParameterRange(m_scalar, 0.0, 10.0);

    m_vector = addVectorParameter("vector", "a vector parameter", ParamVector(1.0, 2.0, 3.0));

    m_choice = addIntParameter("choice", "choose one of the options", Option1, Parameter::Choice);
    V_ENUM_SET_CHOICES(m_choice, Option);
}

MyModule::~MyModule()
{}

bool MyModule::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
    // ports
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

    // parameters
    bool boolean = m_boolean->getValue();
    Float scalar = m_scalar->getValue();
    ParamVector vector = m_vector->getValue();
    Option choice = (Option)m_choice->getValue();

    return true;
}
```
</details>

## How to write a reader module
- create a MyReader test module similar to MyModule above...
- this is the bare minimum also talk about other classes in vistle::Reader (examine, prepareRead, finishRead)

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

`MyModule/CMakeLists.txt`:
```cmake
add_module(MyModule "Test module" MyModule.h MyModule.cpp)
```


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
