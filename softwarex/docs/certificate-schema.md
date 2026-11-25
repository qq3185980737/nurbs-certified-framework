# Certificate schema

A certificate is a UTF-8 JSON object. Numeric geometry is encoded as exact rational
strings:

- `"n/d"`, or
- hexadecimal floating-point strings such as `"0x1.0000000000000p-1"`, which are
  converted exactly with `Fraction.from_float(float.fromhex(...))`.

The verifier never uses binary floating-point for a geometric decision.

## Top level

```json
{
  "tau": "0x1.0000000000000p-2",
  "patches": [ { "nu": 2, "nv": 2, "poles": [ ... ] }, { ... } ],
  "tree": { "box": [ ... 8 values ... ], "kind": "node", "axis": 0, "children": [ ... ] },
  "certified": [ ... ],
  "claims": [ ... ]
}
```

Only `patches` and `tree` are required for the basic coverage certificate. Other
blocks are optional and are checked when present.

## `patches`

Two rational Bézier patches, each with:

```json
{
  "nu": 2,
  "nv": 2,
  "poles": [
    [ "wx", "wy", "wz", "w" ],
    ...
  ]
}
```

`poles` is a row-major `(nu+1) x (nv+1)` control net. Coordinates are homogeneous:
`(wx, wy, wz, w)`. All weights must be strictly positive.

## `tree`

A binary subdivision tree over the 4D parameter box `[0,1]^4`.

```json
{
  "box": [ "u_lo", "u_hi", "v_lo", "v_hi", "s_lo", "s_hi", "t_lo", "t_hi" ],
  "kind": "node",
  "axis": 0,
  "children": [ { ... }, { ... } ]
}
```

A leaf uses:

```json
{ "box": [ ... ], "kind": "leaf", "status": "excluded" }
```

or

```json
{ "box": [ ... ], "kind": "leaf", "status": "output" }
```

The verifier checks that each child box is the exact dyadic split of its parent at the
midpoint of `axis`, and that all boxes lie in `[0,1]^4`.

## `certified`

Optional Krawczyk slices. Each slice is checked independently with exact rational
interval arithmetic. A slice is accepted only if the verifier reproduces the enclosure
and the Krawczyk inclusion test succeeds.

## `claims`

Optional claims for `coverage`, `miss`, `branch_coverage`, and related rules. Claims are
checked only when present; an absent claim is not treated as proved.

## Chunked certificates

A certificate may contain a `chunks` array instead of a single tree. Use:

```bash
ssi-cert verify cert.json --chunked
```

The verifier checks the chunk partition, each chunk's frontier, and its content before
accepting the combined certificate.

## What `ACCEPT` means

`ACCEPT` means the certificate satisfies the implemented contract. It does not mean the
kernel that generated the geometry is correct, nor that the geometric statement is true
outside the stated contract.
