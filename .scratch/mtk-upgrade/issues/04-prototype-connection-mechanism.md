# Prototype the bond connection mechanism on MTK 11

Type: prototype
Status: claimed
Blocked by: 03

## Question

Ticket 03 recommends the hybrid design (0-junction = native `connect` on a `Power` port with `Flow` flow; 1-junction = component with ±1-signed ports; TF/GY as two-port components) and already ran RC, DC motor and 2-DoF against the golden data. Does it hold for the remaining examples (N-DoF, Duffing, a modulated two-port, a custom nonlinear `@oneport`-style element), and does model-building code written with it keep today's ergonomics? A throwaway prototype (starting from the scripts referenced in research/03) should express all six examples in the proposed user-facing API and reproduce `test/regression/reference/`, so the maintainer can react to how model-building code looks with it.
