# RaiSin file formats (v1.3)

All RaiSin data files are UTF-8-compatible, line-oriented text. Numeric fields
are separated by spaces or tabs. Vertex and part identifiers are zero-based in
the CLI formats described here.

## Directed red-black hypergraph (`.rzn2`)

The first line contains five integers:

```text
<pins> <vertices> <hyperedges> <red_vertices> <weight_dimensions>
```

It is followed by exactly `<hyperedges>` lines. Each line represents one
directed hyperedge:

```text
<hyperedge_weight> <source_vertex> <sink_vertex> [<sink_vertex> ...]
```

The source is the first vertex identifier after the hyperedge weight. Repeated
vertex identifiers inside one hyperedge are ignored, but still count toward the
declared raw pin count. All identifiers must be in `0..<vertices>`.

The final `<vertices>` lines describe vertices in identifier order:

```text
<red> <delay> <right_criticality> <weight_0> ... <weight_n>
```

- `<red>` is `0` for combinational logic or `1` for a register/external I/O;
- `<delay>` is a non-negative integral value; forms such as `552.0` are accepted;
- `<right_criticality>` is a non-negative integer;
- the number of resource weights must equal `<weight_dimensions>`;
- all resource weights are non-negative integers.

The loader rejects missing lines, extra vertex fields, out-of-range identifiers,
inconsistent pin/red counts and invalid numeric fields. After parsing, the CLI
also validates red masks, incidence data and combinational acyclicity.

Example:

```text
3 3 1 1 1
10 0 1 2
1 100.0 200 1
0 50.0 100 2
0 20.0 20 1
```

## Target architecture (`.arch`)

The first line contains the number of parts, the number of undirected
connections and, optionally, the number of resource dimensions:

```text
<parts> <connections> [<resource_dimensions>]
```

Each following connection line contains:

```text
<capacity> <delay> <part_u> <part_v>
```

`capacity` is the maximum number of distinct hyperedge signals that may use the
physical link. RaiSin routes every source-to-sink connection over a deterministic
lowest-delay path. A multi-sink hyperedge consumes one unit on every link in the
union of its routes, even when several sinks share a route prefix.

Partition-producing modes attempt deterministic capacity repair before writing
a solution. Communication repair strictly reduces routed overload while
preserving explicit resource capacities, or the existing peak resource load in
legacy balance mode. `eval` never changes its input and rejects violations.

Capacities and delays must be non-negative. Part identifiers must be in
`0..<parts>` and a connection cannot link a part to itself. The loader mirrors
each connection, so the architecture matrix is symmetric. Missing connections
are routed over the lowest-delay path and remain distinct from explicitly
declared zero-delay connections. Architectures must be connected; duplicate
connections and disconnected topologies are rejected.

When `<resource_dimensions>` is present and greater than zero, the connection
lines are followed by exactly one capacity row per part, ordered by part ID:

```text
<part_id> <capacity_0> ... <capacity_n>
```

The number of capacities on every row must equal `<resource_dimensions>`. This
dimension must also match the vertex-weight dimension of the `.rzn2` circuit
when the architecture is used for partitioning or evaluation. Omitting the third
header field preserves the legacy balance-factor behavior.

For an explicit matrix, RaiSin computes each part/resource load as the sum of
the corresponding vertex weights and requires it to be less than or equal to the
declared capacity. This is constraint 5.2f of the thesis.

Example:

```text
3 3
10000 300 0 1
10000 300 0 2
10000 300 1 2
```

## Partition solution (`.sol`)

A solution contains exactly one non-negative part identifier per vertex, in
vertex order:

```text
0
1
1
0
```

The default build stores part identifiers as `uint8_t`, so values are limited to
`0..255`; a given run additionally expects them to be below `part_number`.
RaiSin rejects truncated files, malformed values and identifiers outside the
compiled `PART` range.

The `cluster` mode writes integer cluster identifiers without truncation, because
intermediate cluster counts may exceed the placed-partition limit. Such a file is
loadable by `refine` or `eval` only when every identifier is within the compiled
`PART` range and below the requested `part_number`.

The `partfile` CLI option is an output prefix, not the final name: the writer
appends `.sol`. In `refine` and `eval`, `part_file`/`partfile` name an existing
solution and are read literally as supplied.
