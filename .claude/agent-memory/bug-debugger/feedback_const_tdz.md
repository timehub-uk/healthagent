---
name: const-tdz-before-use
description: `const`/`let` declarations inside forEach/loop callbacks used before they are declared in the same callback scope — temporal dead zone crash
type: feedback
---

In `healthagent/ui/templates/index.html`, the `rungs.forEach` callback in `drawHelix` declared `isPurine1` with `const` on a line AFTER it was referenced. JavaScript `const` (and `let`) are block-scoped and NOT hoisted to usable values — accessing them before the declaration line throws `ReferenceError: Cannot access 'isPurine1' before initialization`.

**Why:** The variable was inside a `const` declaration lower in the same callback body. Unlike `var` or `function` declarations, `const`/`let` sit in the Temporal Dead Zone from the top of their block until their declaration line is reached.

**How to apply:** When reviewing JS in this codebase, scan for any `const`/`let` variable used before its declaration line within the same block scope, especially inside array callbacks (`forEach`, `map`, etc.) where the ordering can be easy to miss.
