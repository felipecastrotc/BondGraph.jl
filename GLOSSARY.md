# BondGraph.jl

A bond-graph dialect of ModelingToolkit: bond-graph models are ordinary MTK systems, so they compose, couple, and simulate with the same tools as any other MTK model.

## Language

### Power and bonds

**Power port**:
The connector through which an element exchanges power, carrying one effort and one flow.
_Avoid_: pin, connector (when meaning this specifically)

**Effort**:
The intensive power variable of a port (force, voltage, pressure).
_Avoid_: potential, across variable

**Flow**:
The extensive power variable of a port (velocity, current, volume flow).
_Avoid_: through variable

**Bond**:
A connection between two power ports, with a direction (the half-arrow) that fixes the positive sense of power transfer. An element's law is written in the variables of its bond as drawn.
_Avoid_: link, edge, wire

### Elements

**Element**:
A system with one or more power ports whose constitutive law relates their efforts and flows.
_Avoid_: component, block

**One-port element**:
An element with a single power port: a source, an inertia, a compliance, or a resistance.

**Source**:
A one-port element that imposes an effort (Se) or a flow (Sf).

**Inertia (I)**:
A one-port element that stores energy in its flow's conjugate momentum; mechanically a mass, electrically an inductance.
_Avoid_: mass (outside mechanical examples)

**Compliance (C)**:
A one-port element that stores energy in its displacement; mechanically a spring, electrically a capacitance.
_Avoid_: spring, capacitor (outside domain-specific examples)

**Resistance (R)**:
A one-port element that dissipates power; mechanically a damper, electrically a resistor.
_Avoid_: damper (outside mechanical examples)

**Custom element**:
A user-defined element, often nonlinear, built with the same ergonomics as the built-in ones.

### Junctions and two-ports

**0-junction**:
A node where all bonds share one effort and their flows sum to zero, each flow counted positive when its bond points into the junction.

**1-junction**:
A node where all bonds share one flow and their efforts sum to zero, each effort counted positive when its bond points into the junction.

**Transformer (TF)**:
A two-port element that scales effort and flow by reciprocal factors, conserving power.

**Gyrator (GY)**:
A two-port element that relates the effort of each port to the flow of the other, conserving power.

**Modulated two-port (mTF, mGY)**:
A transformer or gyrator whose ratio is an expression of system variables rather than a constant.
