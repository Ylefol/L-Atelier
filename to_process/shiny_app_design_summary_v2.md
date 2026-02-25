# Shiny Workflow Builder — Design Summary (v2)

## Context

This document summarizes design decisions for an R Shiny application that serves as a guided front-end for an existing R toolkit for bulk sequencing analysis. The toolkit is largely complete and covers functionality such as differential gene expression, ATAC-seq peak analysis, QC, and normalization. The Shiny app is a new layer being built on top of it.

**Target audience:** Wet lab biologists with little to no bioinformatics or R experience. The app is not intended to replace the need for developing those skills, but to scaffold their development — particularly through its script export feature (see below).

**Implementation constraint:** The app should be buildable in pure R and Shiny with no custom JavaScript. This is a deliberate decision to keep the codebase maintainable and debuggable by someone without JavaScript experience.

---

## Core Concept

The app is a **constrained, guided workflow builder**. Users construct an analysis pipeline by placing functional "boxes" (each representing a function from the toolkit) and connecting them. The key design principle is that the app **restricts available options based on context** — if a user places an "Input RNA counts" box, only functions that accept raw RNA counts as input are offered as the next step. This prevents incompatible chaining and reduces decision paralysis for non-expert users.

---

## Workflow Logic

- The underlying model is a **directed acyclic graph (DAG)**: nodes are toolkit functions, edges represent compatible input/output relationships.
- **Branching is supported**: a single box can have multiple boxes bound to its output. For example, an "Input RNAseq" box could have both a QC box and a Normalization box connected to it, running in parallel.
- When a user adds a box, the app presents a **parameter dialog** for that function, exposing only the parameters relevant to the user (things like output paths and internal plumbing are handled automatically by the app).
- Users must also be able to **revisit and edit parameters of previously placed boxes**.
- Execution is **deferred** — no R code runs during the interactive building phase. Everything executes only when the user finalizes the workflow.

---

## Interaction Model

Adding a new step follows a three-stage modal dialog:

1. **What step?** — the user selects a function from a list filtered to only those compatible with what already exists in the workflow.
2. **Connect to what?** — the user selects which existing box feeds into the new one. This implicitly determines the new node's position in the graph, so no drag-and-drop or manual positioning is needed.
3. **Parameters?** — the user tunes the relevant parameters for that function (see config-driven design below).

This approach avoids any need for a drag-and-drop canvas while remaining intuitive, since it mirrors how a workflow actually works — the user is deciding what feeds into what, which is the meaningful decision.

---

## Visual Design

- The workflow is rendered as an **interactive node graph using the `visNetwork` R package**, which is natively compatible with Shiny and requires no JavaScript.
- The graph uses **`visHierarchicalLayout`** to produce a structured top-down layout (similar in appearance to the Galaxy workflow editor), rather than the default physics-based floating layout.
- The graph serves primarily as a **visualization of the workflow state** rather than the primary interaction surface. Users build the workflow through the sidebar/modal system described above; the graph updates reactively to reflect what has been built.
- Nodes represent functions, directed edges represent data flow between them.

---

## Script Export

A central feature: once a workflow is complete, the app **writes it out as a clean R script** that the user can download, run independently, and modify. This serves both practical and pedagogical purposes — users can see exactly what their clicks produced in code, supporting the gradual development of R literacy.

---

## Backend Architecture

### Config-Driven Design
Each function in the toolkit should be described by a **structured metadata config** (e.g. YAML or JSON), stored alongside the R package. Each entry should define at minimum:

- Accepted input type(s)
- Output type(s)
- Tunable parameters (name, type, default value, optional description)
- Parameters to hide from the user (e.g. output path, internal flags)

The Shiny app reads from this config dynamically to populate menus, parameter dialogs, and valid next-step options. **No function information should be hardcoded in the Shiny UI logic itself.**

### Updatability
This config-driven approach is specifically chosen to make the app maintainable as the toolkit evolves. Adding a new function, changing a parameter, or modifying input/output types should require only a config update, not changes to Shiny code.

### Execution Model
The interactive session is kept lightweight — the app manages state (which boxes exist, how they're connected, what parameters are set) but does not execute toolkit functions during this phase. Execution happens at the end, triggered explicitly by the user.

---

## Out of Scope (for now)

- Deployment and hosting decisions are pending.
- Drag-and-drop canvas interaction is explicitly **not** being pursued, in favor of the modal-based connection approach described above.

---

## Summary of Design Priorities

1. Complexity lives on the developer's side, not the user's.
2. The app guides rather than replaces expertise.
3. No custom JavaScript — everything is implemented in R and Shiny.
4. The config is the single source of truth for what the app knows about the toolkit.
5. The exported R script is a first-class output, not an afterthought.
