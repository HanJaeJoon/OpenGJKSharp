---
title: "Local machine now has .NET SDK 9.0.200 and 10.0.402, so all target frameworks build…"
author: Jaejoon
date: 2026-11-01
id: l-19b3a1b4
scope: repo
status: active
touches:
  - tests/OpenGJKSharp.Tests
---

## Decided
As of 2026-11-01 the local machine has .NET SDK 9.0.200 and 10.0.402 installed, so plain
`dotnet build` / `dotnet test tests/OpenGJKSharp.Tests` work for every target framework
(netstandard2.0;netstandard2.1;net8.0;net10.0 for the library, net8.0;net10.0;net48 for tests).

**Why:** before 2026-11-01 only SDK 9.0.200 was present and net10.0 failed with NETSDK1045, so
earlier notes and dogfooding rounds restricted builds with `-p:TargetFrameworks=net8.0`.

## Watch out
no TFM override is needed any more; drop `-p:TargetFrameworks=...` from commands
copied out of older notes. Re-check `dotnet --list-sdks` if a build fails with NETSDK1045.

_From Claude Code's memory of this repo, last updated Nov 1._
