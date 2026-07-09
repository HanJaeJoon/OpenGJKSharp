## OpenGJKSharp
OpenGJKSharp is a C# native implementation of the GJK (Gilbert-Johnson-Keerthi) algorithm, designed for efficient collision detection between convex polyhedra in 3D space.
This project is inspired by the original [openGJK](https://github.com/MattiaMontanari/openGJK) (written in C) and reimagined for the .NET ecosystem.

## Installation
You can install `OpenGJKSharp` via NuGet:

``` bash
dotnet add package OpenGJKSharp --version 0.3.0
```

Or using the Package Manager:

``` bash
Install-Package OpenGJKSharp -Version 0.3.0
```

> **Note:** Starting with 0.3.0, the entry-point class is named `OpenGjk`
> (previously `OpenGJKSharp`, which collided with the namespace and could not be
> referenced without a namespace alias).

## Usage

Here is an example of detecting a collision between two overlapping cubes:

``` c#
using System.Numerics;
using OpenGJKSharp;

// Cube 1
var a = new Vector3[]
{
    new(0, 0, 0),
    new(1, 0, 0),
    new(0, 1, 0),
    new(1, 1, 0),
    new(0, 0, 1),
    new(1, 0, 1),
    new(0, 1, 1),
    new(1, 1, 1),
};

// Cube 2
var b = new Vector3[]
{
    new(0.5f, 0.5f, 0),
    new(1.5f, 0.5f, 0),
    new(0.5f, 1.5f, 0),
    new(1.5f, 1.5f, 0),
    new(0.5f, 0.5f, 1),
    new(1.5f, 0.5f, 1),
    new(0.5f, 1.5f, 1),
    new(1.5f, 1.5f, 1),
};

bool hasCollision = OpenGjk.HasCollision(a, b);

Console.WriteLine($"Collision detected: {hasCollision}"); // Outputs: true
```

The following image illustrates the collision between the two cubes:

![Collision Example](https://raw.githubusercontent.com/HanJaeJoon/OpenGJKSharp/refs/heads/main/example.png)

### Penetration depth and contact normal (EPA)

`ComputeCollisionInformation` runs the GJK algorithm and, when a collision is detected,
the EPA (Expanding Polytope Algorithm) to compute the penetration depth and contact normal:

``` c#
double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 contactNormal);

if (distance < 0)
{
    // Colliding: -distance is the penetration depth, contactNormal is the contact normal
    Console.WriteLine($"Penetration depth: {-distance}, normal: {contactNormal}");
}
else
{
    // Separated: distance is the minimum distance between the two bodies
    Console.WriteLine($"Distance: {distance}");
}
```

### Minimum distance and closest points

`ComputeMinimumDistance` returns the minimum distance between two bodies
(0 when they collide) and can also report the closest point on each body:

``` c#
double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

Console.WriteLine($"Distance: {distance}");
Console.WriteLine($"Closest points: {closestA} <-> {closestB}");
```

`HasCollision` and `ComputeMinimumDistance` also accept `Vector2[]` polygons
(treated as flat shapes on the z = 0 plane).