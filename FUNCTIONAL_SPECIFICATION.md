# ![#f03c15](https://placehold.co/15x15/E04606/E04606.png) Functional specification

The functional specification for the Polytunnel-PV module.

#### Table Of Contents

1. [Requirements](#1-requirements)
    1. [General requirements](#11-general-requirements)
    2. [Interface requirements](#12-interface-requirements)
    3. [Logging requirements](#13-logging-requirements)
2. [Use cases](#2-use-cases)
    1. [Simulating performance](#21-simulating-performance)
    2. [Optimising performance](#22-optimising-performance)
        1. [Design parameters](#221-design-parameters)
        2. [Operation parameters](#222-operation-parameters)

# 1 Requirements

Outlined are the requirements of the Python module.

# 1.1 General requirements

The code should be able to simulate curved photovoltaic (PV) modules over a range of curved topologies with different PV-cell current-voltage (IV) curve characteristics,
with and without bypass diodes present,
and give an assessment of the performance of each individual cell, as well as the overall module performance,
both when operating at maximum power point (MPP) or when operating under different conditions.

# 1.2 Interface requirements

There should be multiple interfaces through which users and developers can engage with the module:

* A command-line interface (CLI) for the running of the code, simulating the performance of curved modules etc.;
* Exposed APIs for determining the performance of such modules for running within [CLOVER](https://github.com/CLOVER-energy/CLOVER)
  * **NOTE,** this will mean that code from CLOVER cannot be imported in this package to avoid circular dependencies;
* (Ideally) integration with the [CLOVER GUI](https://github.com/CLOVER-energy/CLOVER-GUI).

# 1.3 Logging requirements

Good logging at debug, info, warning, error and critical levels throughout.

# 2 Use cases

There are two main use cases, whether directly accessed or through an API: simulation and optimisation.

## 2.1 Simulating performance

Given some orientation, configuration of bypass diodes, cell geometry, module geometry, etc., a user could likely wish to consider how such a PV module may perform under some environmental conditions (ambient temperature, wind speed, irradiance etc.).

## 2.2 Optimising performance

A user may also wish to optimise the performance of such modules. This may be through optimisation at the design or operation stage. **Note,** these optimisations may either take place within the Python module or externally, utilising its functionality and being supported by it. This will be determined at the design stage.

## 2.2.1 Design parameters

A user may wish to optimise any number of parameters:

* The number of cells in the module, given some length of the module;
* The number of, and position of, bypass diodes within the module;
* The cut-off voltage of the bypass diodes included in the module.

## 2.2.2 Operation parameters

Once the above design parameters have been fixed, either as the result of a user-instigated optimisation or through a design consideration, a user may wish to optimise any number of operation parameters:

* The orientation of the module (and, by extension, of the polytunnel);
* The direction of curvature, i.e., whether the cells are arranged along the length of the polytunnel (each experiencing the same irradiance) or around an arc of the polytunnel (with each receiving a different irradiance).
