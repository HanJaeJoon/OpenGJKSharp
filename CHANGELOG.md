# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `netstandard2.1` target: the package can now be used from .NET Core 3.0+, Unity 2021.2+, and Mono in addition to .NET 8
- SourceLink and symbol package (`.snupkg`) for source-level debugging

### Changed

- NuGet publishing now happens only on `v*` release tag push (previously on every push to `main`), with a tag-vs-csproj version check

## [0.3.0] - 2026-07-10

### Changed

- **Breaking:** rename the entry-point class `OpenGJKSharp` to `OpenGjk`. The old name collided with the namespace, so consumers could not reference it without a namespace alias
- Replace `Console.WriteLine` noise in the library with `Debug.WriteLine`
- Ship IntelliSense XML documentation with the package

### Added

- Input validation with clear exceptions: null/empty vertex arrays, non-finite (NaN/Infinity) coordinates, negative precision, and `GkPolytope` shape mismatches
- `ComputeMinimumDistance` overloads for `Vector3[]`/`Vector2[]` that return the distance and the closest point on each body (witness points)

### Fixed

- EPA reporting zero penetration depth when both bodies share their first vertex
- EPA mis-orienting face normals when a polytope face passes through the origin, which stalled EPA at ~0 penetration depth for sphere-like hulls
- Missing returns in degenerate witness branches that could produce NaN results

## [0.2.0] - 2026-07-09

### Changed

- Apply modern C# conventions to the model types (public fields to properties)
- Fix the `_presicion` typo to `_precision`

## [0.1.0] - 2026-07-09

### Added

- EPA (Expanding Polytope Algorithm) ported from upstream openGJK: `ComputeCollisionInformation` returns the penetration depth and contact normal
- `Vector3` convenience overload and a `ComputeMinimumDistance` overload accepting a `GkSimplex`

### Changed

- NuGet publishing switched from an API key to Trusted Publishing (OIDC)

## [0.0.5] - 2025-03-08

### Added

- Initial release: GJK minimum distance and `HasCollision` for convex polyhedra (`Vector3[]`) and flat polygons (`Vector2[]`)

[Unreleased]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.0.5...v0.1.0
[0.0.5]: https://github.com/HanJaeJoon/OpenGJKSharp/releases/tag/v0.0.5
