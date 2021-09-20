Organization of the code
========================

The MOM6 source code is divided into a tree of directories (folders) to group
related code (e.g. `src/core`) or similar modules (e.g.
`src/parametizations/vertical`).

The highest level of directories are:

`src/`
  Code underneath `src/` is always required and compiled.

`config_src/`
  Under `config_src` are various drivers and memory configuration sources that
  can only be compiled in limited configurations. See :ref:`config_src`

`pkg/`
  Packages (source code) from other sources/parties only some of which might
  be used. We include the entire package as a sub-module but use
  symbolic-links to extract the parts the MOM6 uses.

`docs/`
  The directory that contains this documentation, namely that beyond the
  in-code API documentation. Some of the files are configuration files
  needed for running doxygen and sphinx. Most documentation in this folder
  is in the form of reStructuredText (.rst) files.

The directory tree is::

  MOM6
  ├── config_src
  │   ├── coupled_driver
  │   ├── dynamic
  │   ├── dynamic_symmetric
  │   ├── ice_solo_driver
  │   ├── solo_driver
  │   └── unit_drivers
  ├── docs
  │   └── images
  ├── pkg
  │   ├── CVMix-src
  │   │   ├── ...
  │   │   └── src
  │   │       ├── drivers
  │   │       └── shared
  │   └── GSW-Fortran
  └── src
      ├── ALE
      ├── core
      ├── diagnostics
      ├── equation_of_state
      │   └── TEOS10
      ├── framework
      ├── ice_shelf
      ├── initialization
      ├── parameterizations
      │   ├── CVmix -> ../../pkg/CVMix-src/src/shared
      │   ├── lateral
      │   └── vertical
      ├── tracer
      └── user

.. _config_src:

`config_src/`
-------------

`dynamic/`, `dynamic_symmetric/`
  One or none of `config_src/dynamic/` or `config_src/dynamic_symmetric/` can
  be included at compile time. If neither is used then a `MOM_memory.h` file
  specific to the model configuration must be present - this is known as a
  "static" compile with fixed layout and domain shape.
  
`solo_driver/`
  This driver produces an ocean-only executable with no other coupled
  components (no sea-ice, no atmosphere, etc.). It is the simplest
  configuration and fastest to compile and thus used for a lot of testing.

`coupled_driver/`
  This driver provides an interface for the GFDL coupler to call. It requires
  compiling MOM6 along with at least a sea-ice model and possibly all other
  components in a coupled model.

.. _src:

`src/`
------

`core/`
  The dynamical core modules (except for the ALE remapping/regridding).

`ALE/`
  Functions for remapping from between arbitrary vertical grids
  and generating grids.

`diagnostics/`
  Some diagnostic calculations

`equation_of_state/`
  Various equations of state (linear; Wright, 1997; TEOS-10; ...).

`framework/`
  Low-level wrappers for communication, diagnostics management, parsing
  of input parameters, time management, CPU clocks.

`initialization/`
  Initialization of the horizontal grid, vertical coordinate, and the state.

`parameterizations/lateral`
  Sub-grid scale parameterization with fluxes primarily oriented in the
  lateral direction.

`parameterizations/vertical`
  Sub-grid scale parameterization with fluxes primarily oriented in the
  vertical direction, including the top and bottom boundary layer schemes.

`tracer/`
  Everything to do with tracers, including advection and isopycnal stirring.

`user/`
  Initialization and forcing for specific (coded) configurations.
