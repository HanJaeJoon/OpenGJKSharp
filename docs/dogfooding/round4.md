# OpenGJKSharp 4차 dogfooding (0.7.0-rc.1)

- 날짜: 2026-11-01
- 대상 패키지: `OpenGJKSharp 0.7.0-rc.1` (`main` 2e4ebb7에서 `dotnet build -c Release -p:Version=0.7.0-rc.1`로 만든 로컬 피드, TFM: netstandard2.0 / netstandard2.1 / net8.0 / net10.0)
- 소비자 프로젝트 (저장소 밖 `../OpenGJKSharp-dogfood4/`, 라운드 종료 후 삭제)
  - `Dogfood4/` (net8.0) - 시나리오 1~4, 7. `dotnet run -c Release -- <1|2|2c|3|3p|4|7>`
  - `Dogfood4/Oracle.cs` - 라이브러리와 독립적인 double 해석 오라클 (구 이차식, 박스 slab, 캡슐 선분-구, 볼록 다면체 평면 클리핑, 정확한 박스 SDF, 2D SAT)
  - `Upgrade060/v060`, `Upgrade060/v070rc1` (net8.0) - 같은 `Program.cs`를 nuget.org의 0.6.0과 0.7.0-rc.1에 빌드. 시나리오 4의 시간 측정과 시나리오 5
  - `NetFxConsumer/` (net48) - 시나리오 6
  - net10.0은 `Dogfood4`를 `-p:TargetFramework=net10.0`으로 다시 빌드해 시나리오 7과 P1 재현을 돌렸다
- 이번 라운드 초점: 0.7.0 신규 API (`pointA`/`pointB` 접촉점, `TryCastRay`/`TryCastShape`, `HasCollision` 조기 종료)
- 시나리오는 패키지에 들어 있는 README와 XML 문서만 보고 작성했다. 소스는 P1 원인을 진단할 때만 열었다(`OpenGjk.Cast.cs`).

## 결론

**P0 없음, P1 1건.** ray cast와 shape cast의 수렴 허용 오차가 world 좌표 크기에 비례해서, 같은 도형이라도
원점에서 멀어질수록 결과가 나빠진다. 원점 근처에서는 모든 시나리오가 오라클과 float 양자화 수준으로
일치했지만, 같은 박스를 x = 1000에 두면 윗면에서 1/1024 떨어져 지나가는 광선을 hit로 보고하고,
평면 법선이 최대 0.1 (약 6도) 틀어진다. 보고된 fraction이 실제 충돌 시점을 넘는 일은 없어 터널링은 생기지 않는다.
0.7.0 릴리스 전에 고치거나, 최소한 README/XML의 "exact on flat faces", "about 1e-3" 문구에 좌표 크기 조건을 붙일 것을 권고한다.

나머지는 P2 2건(그중 1건은 0.6.0부터 있던 EPA 문제), friction 3건, 미검증 2건이다.
접촉점, `HasCollision` 동치, 0.6.0 -> 0.7.0 출력 비교, net48, net10.0, workspace 할당 0은 모두 통과했다.

## 시나리오별 결과

### 1. ray cast vs 해석 오라클 [원점: 통과 / 원점 밖: P1]

도형 6종(구, AABB, 회전 박스, 캡슐, 정팔면체 `ConvexHull`, 회전 박스 `ConvexHull`)마다 원점 근처 무작위 광선
3000개(10%는 내부 시작), 좌표 (1e4, -1e4, 1e4) 근처 1000개를 쐈다. 짝수 번째는 workspace 없이, 홀수 번째는 workspace로 호출했다.
오차는 거리 단위(fraction 차이 x `|translation|`)다.

원점 근처:

| 도형 | hit/miss 불일치 | hit 지점 오차 max | 실제 hit보다 뒤 max | 법선 오차 max (p50) |
| --- | --- | --- | --- | --- |
| 구 | 0 / 3000 | 2.1e-5 | 8.7e-8 | 1.95e-3 (8.4e-4) |
| AABB | 0 / 3000 | 6.8e-6 | 1.4e-6 | 8.5e-12 |
| 회전 박스 | 0 / 3000 | 6.5e-6 | 8.7e-7 | 3.9e-7 |
| 캡슐 | 0 / 3000 | 1.0e-4 | 1.1e-7 | 4.6e-3 (9.9e-4) |
| 정팔면체 hull | 0 / 3000 | 2.5e-11 | 2.5e-11 | 5.1e-8 |
| 회전 박스 hull | 0 / 3000 | 1.9e-6 | 2.5e-7 | 2.9e-7 |

- "실제 hit보다 뒤" 최대 1.4e-6은 support point의 float 양자화 때문이다(`Box`는 꼭짓점을 float로 계산한다).
  꼭짓점이 float로 정확히 표현되는 박스로 다시 재면 1.6e-15다(아래 P1 표).
- 곡면 법선 오차는 구가 최대 1.95e-3, 캡슐이 최대 4.6e-3 (p50은 둘 다 ~1e-3)이다. 문서의 "about 1e-3"보다 캡슐 최대값이 4.6배 크다.
  원인은 P1과 같은 허용 오차로 보인다(광선 시작점까지 포함한 좌표 크기 ~10에서 허용 오차가 ~1e-5가 된다).
- 특수 케이스는 전부 문서대로 동작했다:
  - README 예제 박스: `0.4`, `<-1, 0, 0>`
  - 구에 접하는 광선(y = 1): hit, t = 0.49990. y = 1 - 1e-4: t = 0.4985854 (해석값 0.4985858). y = 1 + 1e-4: miss
  - 표면에서 시작(안쪽, 바깥쪽, 접선 방향 모두), 내부에서 시작, 내부/표면에서 translation 0: hit, fraction 0, 법선 0
  - 외부에서 translation 0, 반대 방향, 닿기 전에 끝나는 광선: miss. 끝점이 정확히 면 위: hit, t = 1
  - 길이 1e6 광선: t = 4e-6, 반지름 1e-3 구: t = 0.39999998, 1e4에 둔 박스/구 정면: t = 0.4
  - 평면 사각형 `ConvexHull` 앞/뒤: `0.4`, `<-1, 0, 0>` / `0.6`, `<1, 0, 0>`. 선분, 점 hull: `0.4`
  - NaN origin, 무한 translation: `ArgumentException`
  - 예외 하나: 박스 윗면과 평행하게 9.5e-7 위를 지나는 광선이 hit로 나온다. 원점 근처에서도 허용 오차가 ~1e-5라서다(P1)

(1e4, -1e4, 1e4) 근처에서는 도형마다 1000개 중 1~4개가 "라이브러리 hit, 오라클 miss"였고(가장 가까이 지나간 거리 최대 1.3e-2),
평면 법선 오차가 최대 1.4 (회전 박스 1.33, 정팔면체 1.16)까지 나왔다. 구의 법선 오차는 p50 5.1e-2, 최대 0.17이다.
float 양자화(1e4에서 간격 ~1e-3)만으로는 설명되지 않아 P1로 따로 분석했다.

### 2. shape cast vs 오라클 [원점: 통과 / 원점 밖: P1]

쌍 5종: 구-구(해석), 구-회전 박스(정확한 박스 SDF 위에서 이분법), AABB-AABB(Minkowski 박스에 slab),
캡슐-회전 박스와 hull-hull(`ComputeMinimumDistance` 위에서 이분법, 라이브러리 거리 기준).

원점 근처:

| 쌍 | 불일치 | 시작 겹침 | toi 오차 max | 실제보다 뒤 max | 법선 오차 max | toi로 옮긴 뒤 depth max | toi + 1e-3이면 겹침 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 구-구 | 0 / 2000 | 213 | 1.1e-5 | 5.8e-7 | 1.8e-3 | 3.7e-7 | 1390 / 1390 |
| 구-회전 박스 | 0 / 2000 | 306 | 3.3e-5 | 1.3e-6 | 4.6e-3 | 3.4e-7 | 1416 / 1416 |
| AABB-AABB | 0 / 2000 | 341 | 5.8e-7 | 5.8e-7 | 8.1e-12 | 2.4e-7 | 1390 / 1390 |
| 캡슐-회전 박스 | 0 / 400 | 33 | 2.7e-5 | 3.4e-7 | - | 3.1e-7 | - |
| hull-hull | 0 / 400 | 53 | 2.8e-7 | 2.8e-7 | - | 1.9e-7 | - |

- 접촉점 `point`가 b 표면 위에 있는지(`|sdf_b|`): 구-구 최대 4.4e-6, 구-박스 1.4e-6, AABB 2.4e-7
- 시작 겹침은 전부 fraction 0, 법선 0이었고, `point`는 두 도형 모두의 내부였다(문서대로)
- 같은 시나리오를 (1e3, 0, -1e3)으로 옮기면 구-구 법선 오차 p50 1.3e-2 (최대 3.3e-2), 구-박스 toi 오차 최대 1.4e-2,
  AABB-AABB 면 내부 법선 오차 최대 1.0이 나왔다(P1). toi로 옮긴 뒤 1e-5 넘게 겹친 쌍(최대 9.7e-5)도 생겼지만,
  1e3에서 float 간격이 6.1e-5이므로 이쪽은 입력 양자화 범위다.

특수 케이스:

| 구성 | 결과 |
| --- | --- |
| README 예제 (구 -> 박스) | t = 0.3, `<1, 0, 0>`, point `<4, 0.00023, 0.00023>` (README 주석은 `(4, 0, 0)`) |
| 박스-박스 면 대 면 | t = 0.3, `<1, 0, 0>`, point `<4, 1, 1>` (접촉면의 꼭짓점) |
| 접촉한 채 시작: 멀어지는 방향 / 접선 방향 / 들어가는 방향 | 셋 다 hit, t = 0, 법선 0 |
| 간격 1e-4, 1e-5, 1e-6에서 접선 방향 | miss |
| 평행하게 y = 2.001 (miss) / y = 2 (스침) | miss / hit t = 0.39975 |
| translation 0, 떨어짐 / 겹침 | miss / hit t = 0 |

**캐릭터 컨트롤러.** 캡슐(반지름 0.3, 반길이 0.4)이 떨어져 바닥(y = 0)에 서고 (2, 0, 1) 속도로 벽(x = 3)에 비스듬히 부딪혀
z 방향으로 미끄러지는 600프레임. 매 프레임 sweep -> toi까지 이동 -> 남은 이동에서 법선 성분 제거를 최대 4회 반복한다.
"접촉한 채 시작"(법선 0)이 나오면 문서대로 `TryGetPenetration`의 법선을 대신 쓴다.

| 설정 | 최종 위치 (기대 x = 2.7, y = 0.7, z ~ 10) | 접촉한 채 시작한 cast | z가 멈춘 프레임 (200 이후) |
| --- | --- | --- | --- |
| 원점, skin 없음 | (1.013, 0.700, 0.506) | 2279 | 399 |
| 원점, skin 1e-3 | (2.699, 0.700, 9.516) | 0 | 0 |
| (1e3, 0, 1e3), skin 없음 | (1.013, 0.700, 0.506) | 2279 | 399 |
| (1e3, 0, 1e3), skin 1e-3 | (2.698, 0.701, 8.675) | 1139 | 0 |

- 네 경우 모두 프레임 끝에 1e-6 넘게 겹친 적은 없다(fraction이 보수적이라는 문서와 일치)
- skin이 없으면 착지한 순간부터 모든 cast가 "접촉한 채 시작"이 되어 그 자리에 멈춘다(friction)
- skin 1e-3은 원점에서는 충분하지만 좌표 1e3에서는 허용 오차(~1.4e-3)보다 작아 절반 가까운 cast가 다시 접촉 상태가 된다(P1)
- z가 10이 아니라 9.5인 것은 매 스텝 skin만큼 덜 움직이기 때문이다

### 3. 접촉점 (`pointA`/`pointB`) [통과, P2 1건]

도형 쌍 8종 x 겹친 쌍 1500개, `Vector3[]` 1500쌍, `Vector2[]` 2000쌍(SAT 오라클), 1e6 좌표의 `ReadOnlySpan<double>` 1000쌍.

| 검사 | 결과 (전 종류 최대) |
| --- | --- |
| `dot(pA - pB, n) - depth` | 3.2e-7 (span 1e6: 1.6e-10) |
| 각 점이 자기 도형의 support plane 위 (`dot(pA, n) == max dot(A, n)`) | 4.1e-7 (span: 4.7e-10) |
| 각 점이 자기 표면 위 (`|sdf|`) | 2.3e-6 |
| 4인자 `TryGetPenetration`과 depth/법선 동일 | 12000 / 12000 비트 단위 동일 |
| workspace 변형과 동일 | 모든 입력 종류에서 비트 단위 동일 |
| a를 `-(depth + 1e-4) * n`만큼 옮기면 분리 | 100% |
| 2D: depth vs SAT | 9.7e-8, 점이 경계 위 1.4e-7, 옮긴 뒤 SAT depth 2.0e-7 |
| span 1e6 vs `Vector3[]` 원점 | depth 차 3.0e-11 |
| 출력 span 길이 2 | `ArgumentException` "Must hold at least 3 elements (x, y, z), but holds 2." |

- 점은 항상 표면 위에 있고 dot 관계도 맞는다. 다만 면끼리 만날 때 두 점이 "contact region 안에서 마주 본다"는
  문서 문구는 AABB-AABB 2000쌍 중 392쌍에서 맞지 않았다(아래 P2).
- span 점과 `Vector3[]` 점은 최대 1.4 차이가 났지만 둘 다 위 검사를 통과한다. 면 접촉에서는 유효한 점이 여럿이라 그렇다.
- 구-박스 1쌍은 깊이만큼 옮긴 뒤(실제 겹침 4.5e-7) `TryGetPenetration`이 다시 0.077을 보고했다. 4인자 overload와
  0.6.0에서도 같은 값이라 접촉점 기능과는 무관한 기존 EPA 문제다(아래 P2).

**강체 루프.** 회전하는 박스(반치수 0.5 x 0.25 x 0.35)를 바닥에 떨어뜨리고, 프레임마다 접촉점 하나에
법선 impulse(8회 반복, 누적 clamp)와 마찰을 걸고 위치 보정을 했다. 120 Hz, 10초, 마지막 2초를 측정했다.

| 시작 자세 | 접촉점 | 높이 (정지값) | 최대 기울기 | max `|w|` (rad/s) |
| --- | --- | --- | --- | --- |
| 20/10도 기울임 | `pointA` | 0.2468..0.2477 (0.25) | 0.23도 | 0.39 |
| 20/10도 기울임 | `(pA + pB) / 2` | 0.2470..0.2472 | 0.12도 | 0.26 |
| 20/10도 기울임 | `ComputeMinimumDistance` witness (0.6.0까지의 유일한 점) | 0.2474..0.2706 | 2.45도 | 1.24 |
| 평평하게 | `pointA` | 0.2468..0.2476 | 0.22도 | 0.40 |
| 모서리로 45도 | `pointA` | 0.4975..0.4977 (0.5) | 0.06도 | 0.10 |
| 모서리로 45도 | `(pA + pB) / 2` | 0.7049 (모서리로 선 채 유지) | 45도 | 0.03 |

- 새 접촉점은 겹치면 내부로 붕괴하는 witness보다 확실히 낫다(기울기 10배, 각속도 3배 작음). 관통이나 발산은 없다.
- 그러나 완전히 멈추지는 않는다. 바닥면에 놓인 박스의 `pointA`는 프레임마다 바닥면의 네 꼭짓점 중 하나라서
  (`(+-0.5, -0.25, +-0.35)`에 고르게 분포) 단일 접촉점 solver로는 0.2도 안쪽에서 계속 흔들린다(friction).
- 45도 모서리에서 midpoint가 서 있는 것은 완전 대칭인 불안정 평형이 수치적으로 유지된 것이다.

### 4. `HasCollision` 동치와 속도 [통과]

- 8개 overload 모두 기본 precision은 `1E-06`
- 무작위 쌍(도형 6종 조합, 크기 1e-3~1e4, 좌표 0과 1e3), `Vector3[]`, span, `Vector2[]`, 각각 workspace 유무에서
  precision을 {기본값, 0, d, d의 다음 double, 이전 double, d(1 +- 1e-6), d/2, 2d}로 바꿔 `HasCollision == (ComputeMinimumDistance <= precision)`을
  검사했다: **86,794건 중 불일치 0**
- 간격이 기본 precision 근처(`1e-6 x (1 + 0.05k)`, k = -20..20)인 AABB 쌍, 크기 1e-3 / 1 / 1e3: 123건 불일치 0,
  거리도 해석 간격과 같았다

같은 소스를 0.6.0과 0.7.0-rc.1에 빌드한 시간(workspace 재사용, 7회 중 최소, 2회 실행 범위):

| 쌍 | 0.6.0 | 0.7.0-rc.1 |
| --- | --- | --- |
| 떨어진 32점 cloud | 915~940 ns | 447~474 ns |
| 떨어진 구-박스 | 1375~1644 ns | 302 ns |
| 겹친 32점 cloud | 839~857 ns | 847~853 ns |
| 겹친 구-박스 | 413~492 ns | 491 ns |
| 간격 5e-7 구-박스 | 507~627 ns | 660~662 ns |
| (기준) `ComputeMinimumDistance` cloud | 912~916 ns | 924~925 ns |

떨어진 쌍에서 2.0배, 4.5~5.4배 빨라졌다(CHANGELOG는 1.9배, 9.5배. 도형 배치가 달라 비율이 다르다).
겹친 쌍과 간격이 precision 안쪽인 쌍은 조기 종료가 걸리지 않아 같거나 약간 느리다. 문제로 보지 않는다.

### 5. 0.6.0 -> 0.7.0 업그레이드 [통과]

0.6.0에 있던 API만 쓰는 같은 소스를 두 패키지에 빌드해(경고 0개) 출력 1401줄을 왕복 정밀도(`R`)로 diff했다.
**달라진 줄 0개.**

- 도형 400쌍(`Sphere`/`Box`/`Capsule`/`ConvexHull`/`TransformedShape`): 거리, witness, workspace 거리,
  `ComputeCollisionInformation` + 법선, `TryGetPenetration`, `HasCollision` precision 4종
- `Vector3[]` 200쌍과 같은 쌍을 1e6 좌표 span으로, `Vector2[]` 200쌍
- CHANGELOG의 0.7.0 Changed 항목은 `HasCollision`의 속도뿐이고 값은 같다고 하므로 일치한다

### 6. .NET Framework 소비자 [통과]

net48 콘솔 앱(PackageReference)에서 복원/빌드 경고 0개, .NET Framework 4.8.10008.0에서 실행:

- `TryCastRay` (workspace 유무): `True 0.4 <-1, 0, 0>`
- `TryCastShape(ws)`: `True 0.3 <1, 0, 0> <4, 0.0002327955, 0.0002327955>`
- 접촉점, 도형: depth 0.5, `<4.5, ...>` / `<4, ...>`
- 접촉점, `Vector3[]` 사면체: depth 0.17320508, pA `<0.6, 0.2, 0.2>` (면 x + y + z = 1 위), pB `<0.5, 0.1, 0.1>`
- 접촉점, `Vector2[]`: `0.25 <1, 0> <1, 0> <0.75, 0>`
- 접촉점, span 1e6 (`System.Memory` 경로, workspace 유무): depth 0.173205080756888, pA `(1000000.6, 0.2, 0.2)`

### 7. 성능/할당 [통과]

워밍업 2만 회 후 20만 회 측정, net8.0 (net10.0에서도 같은 경향이고 workspace 변형은 0 B):

| 호출 | ns/call | B/call |
| --- | --- | --- |
| `TryCastRay(box, ws)` hit | 1500 | 0 |
| `TryCastRay(sphere, ws)` hit | 1567 | 0 |
| `TryCastRay(hull 32점, ws)` miss | 202 | 0 |
| `TryCastShape(capsule -> box, ws)` | 4219 | 0 |
| `TryCastShape(capsule -> hull 32점, ws)` | 4532 | 0 |
| `TryGetPenetration` + 점 (capsule, box, ws) | 4608 (4인자: 4496) | 0 |
| `TryGetPenetration` + 점 (`Vector3[]` 32점, ws) | 5259 (4인자: 5065) | 0 |
| `TryGetPenetration` + 점 (`Vector2[]`, ws) | 1130 | 0 |
| `TryGetPenetration` + 점 (span 32점, 1e6, ws) | 5717 | 0 |
| `TryCastRay(box)` workspace 없음 | 1471 | 2784 |
| `TryGetPenetration` + 점 (capsule, box) workspace 없음 | 15150 | 109384 |

workspace 변형은 전부 0 B다. 접촉점 추가 비용은 4인자 대비 2~4%다.

## 발견 사항

### P1. cast의 허용 오차가 world 좌표 크기에 비례해서, 원점에서 멀수록 hit/법선이 틀어진다

```c#
// 모든 좌표가 float로 정확히 표현되는 값이라 입력 양자화는 없다.
var box = new Box(new Vector3(1000, 0, 0), new Vector3(1.5f, 0.75f, 1.25f));   // 윗면 y = 0.75

// 윗면과 평행하게 1/1024 위를 지나는 광선: miss여야 한다
OpenGjk.TryCastRay(box, new Vector3(995, 0.75f + 1f / 1024, 0), new Vector3(10, 0, 0),
    out double t, out Vector3 n);
// -> true, t = 0.35001, n = <-0.989, 0.149, 0>   (박스를 x = 0에 두면 false)
```

같은 현상의 다른 형태:

- 박스를 원점에 두고 광선을 x = -1000에서 시작(길이 2000, 윗면 1/1024 위): hit, 법선 `<-1, 0.00075, 0>`. 도형이 아니라 광선 시작점이 멀어도 생긴다
- `TryCastShape`: 박스(반치수 0.5)를 x = 1000 근처에서 다른 박스 윗면 1/1024 위로 평행 이동: hit, 법선 `<0.958, -0.288, 0>`. x = 0에서는 miss

꼭짓점이 정확히 표현되는 AABB(반치수 1.5 / 0.75 / 1.25)를 좌표만 옮겨 같은 상대 광선 2485개(모두 hit)를 쐈을 때:

| 박스 위치 | 보고된 hit 지점과 박스 사이 거리 max | 면 내부 법선 오차 max | hit 지점 오차 max | 실제보다 뒤 max |
| --- | --- | --- | --- | --- |
| 0 | 1.6e-12 | 7.5e-13 | 2.7e-12 | 1.6e-15 |
| 1e2 | 9.9e-5 | 6.9e-4 | 1.6e-4 | 1.6e-15 |
| 1e3 | 1.7e-3 | 1.1e-1 | 7.2e-3 | 1.6e-15 |
| 1e4 | 1.7e-2 | 1.37 | 6.9e-1 | 1.5e-15 |
| 1e5 | 1.7e-1 | 1.37 | 2.85 | 1.4e-15 |

**원인** (진단을 위해 소스를 읽은 결과): `OpenGjk.Cast.cs`의 `Cast`는 `Norm2(v) <= _castTolerance^2 * scale2`일 때 멈춘다.
`_castTolerance = 1e-6`이고 `scale2`는 지금까지 본 support point의 world 좌표 크기 제곱의 최대값이다
(ray cast는 광선 시작점도 들어간다). 그래서 "닿았다"고 판단하는 거리가 `1e-6 x |world 좌표|`가 된다.
주석에는 support point의 float 양자화(상대 6e-8)보다 위에 두려는 의도가 적혀 있지만,

1. 입력이 정확히 표현되는 좌표에서도 똑같이 느슨해지고(위 표),
2. 그 거리 안에서는 모서리 근처를 지나는 광선이 다른 면의 분리 평면에서 멈춰 법선이 90도 넘게 틀어지고,
3. 곡면 법선 오차("about 1e-3")도 이 허용 오차에서 나온다(`sqrt(tol / r)` 수준으로 추정. 1e4에서 구 법선 오차 p50 5e-2).

**영향.** 좌표 100~1000 단위의 장면은 흔하다. 1e3에서 평면 법선 6도, 1.7e-3 떨어진 광선의 false hit,
캐릭터 컨트롤러가 1.4e-3보다 큰 skin을 요구하는 문제가 생긴다(시나리오 2). cast에는 double 입력 경로가 없어서
큰 좌표를 쓰는 사용자가 피해 갈 방법도 없다. fraction은 보수적으로 유지되므로 터널링이나 겹침은 생기지 않는다.

**권고.** 허용 오차를 world 좌표가 아니라 문제의 국소 크기(광선 점 기준 상대 좌표나 도형 extent)에 맞추거나,
cast 시작 시 기준점(광선 origin이나 b의 center)을 빼고 국소 좌표에서 푼다. 고치지 않고 내보낸다면 README와
XML remarks의 "exact on flat faces", "about 1e-3 for unit-sized shapes"에 "원점 근처에서"라는 조건과
오차가 좌표 크기에 비례한다는 설명을 넣어야 한다.

### P2. 면 접촉의 `pointA`/`pointB`가 "contact region 안에서 마주 본다"는 문서와 다를 때가 있다

```c#
var a = new Box(new Vector3(0.32553646f, 0.15131263f, -0.650346f), new Vector3(0.39351958f, 0.92698365f, 1.7169726f));
var b = new Box(new Vector3(0.42886457f, 1.4347839f, -0.2184791f), new Vector3(0.36379245f, 1.6231194f, 1.3787163f));
OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 n, out Vector3 pA, out Vector3 pB);
// depth 0.65398, n <1, 0, 0>, pA <0.719, 1.078, -0.381>, pB <0.065, 1.180, -0.435>
// 두 박스의 y 겹침 구간은 [-0.188, 1.078]인데 pB.y = 1.180: contact region 밖이고, pA와도 마주 보지 않는다
```

무작위 AABB-AABB 겹친 쌍 2000개 중 1608쌍은 문서대로였고, 387쌍은 두 점이 region 안에 있지만 마주 보지 않았고
(접선 방향 어긋남 최대 0.72), 5쌍은 점 하나가 region 밖이었다. 대부분 깊게 겹친 쌍이다.
표면 위, support plane 위, `dot(pA - pB, n) == depth`는 모든 경우에 성립하므로 depenetration에는 지장이 없다.
README/XML의 "Where flat faces meet, the two points face each other inside the contact region"을
"each on its own surface; on a face contact typically a vertex of the contact region" 정도로 약하게 고치거나,
면 클리핑으로 점을 고르는 것을 권고한다.

### P2 (0.6.0부터). 거의 맞닿은 곡면 쌍에서 EPA 깊이가 크게 과대평가된다

```c#
var s = new Sphere(new Vector3(-1.0468123f, 1.6842266f, 0.11092628f), 1f);
var b = new Box(Vector3.Zero, new Vector3(1.6f, 0.3f, 0.4f),
    new Quaternion(0.21695663f, 0.23592864f, -0.17743732f, 0.9304749f));
OpenGjk.TryGetPenetration(s, b, out double depth, out Vector3 n);
// depth 0.00929, 정확한 박스 SDF로 잰 실제 겹침은 1.75e-6 (5300배)
```

- 0.6.0에서도 같은 값이 나온다(시나리오 5의 diff에 포함). 0.7.0 회귀가 아니다.
- 맞닿게 만든 구-박스 포즈 24,000개 중 1개에서 1e-3 넘게 과대평가했다. 시나리오 3에서는 depenetration 직후 구-박스
  1500쌍 중 1쌍이 0.077을 보고했다(실제 겹침 4.5e-7).
- 문서의 "within about 0.3% for spheres"는 깊이에 대한 상대 오차라 깊이가 0에 가까울 때는 아무것도 보장하지 않는다.
  새 `pointA`/`pointB`도 같은 depth를 쓰므로, 한 번에 depenetration하는 사용자는 위 두 경우에 0.009, 0.077만큼 더 밀게 된다.

### friction. 접촉한 채 시작한 cast는 멀어지는 방향이어도 막힌다

문서("already overlap (or touch) at the start pose")대로지만, 그 결과 toi까지 옮긴 직후의 다음 cast는 방향과 관계없이
fraction 0, 법선 0이 된다. skin 없이 쓴 slide 루프는 착지 순간 멈췄다(시나리오 2, x = 1.01에서 정지, 접촉 cast 2279회).
README에 skin(toi보다 조금 덜 움직이기)을 쓰는 slide 루프 예시와, skin이 cast 허용 오차보다 커야 한다는 점
(현재는 P1 때문에 좌표 크기에 비례)을 적어 두는 것이 좋다.

### friction. 면 접촉의 접촉점은 꼭짓점이라 단일 접촉점 solver는 완전히 멈추지 않는다

바닥에 놓인 박스의 `pointA`는 바닥면의 네 꼭짓점을 오간다(시나리오 3 강체 루프). 쿼리 하나에 점 하나이니 당연한 결과지만,
README의 "Vector3 contact = (pointA + pointB) / 2" 예시만 보면 그 점 하나로 안정한 휴지 접촉이 된다고 기대하기 쉽다.
"휴지 접촉에는 여러 프레임의 점을 모은 manifold가 필요하다"는 한 줄이 있으면 좋다.

### friction. README `TryCastShape` 예제의 `contactPoint` 주석

`contactPoint = (4, 0, 0)`이라고 적혀 있지만 실제로는 `<4, 0.00023, 0.00023>`이다(문서화된 곡면 지연). 주석을 "about (4, 0, 0)"으로 바꾸면 된다.

### 미검증. NuGet package validation

CI에서만 도는 0.6.0 기준 API 검증은 로컬에서 돌리지 않았다. 대신 0.6.0 API만 쓰는 같은 소스가
두 패키지에서 경고 없이 빌드되고 출력이 같다는 것을 확인했다(시나리오 5).

### 미검증. packages.config 방식 복원

라운드 3과 같다. 로컬에 `nuget.exe`가 없어 확인하지 못했다. PackageReference 방식 net48은 확인했다(시나리오 6).

## 이전 라운드에 없던 것

| 항목 | 이전 라운드 | 라운드 4 |
| --- | --- | --- |
| ray cast / shape cast | 없음 (0.7.0 신규) | 광선 24,000개, sweep 13,600개를 독립 오라클과 대조, 캐릭터 컨트롤러 |
| 겹친 상태의 접촉점 | witness만 있었고 겹치면 내부로 붕괴 | 16,500쌍 표면/dot/분리 검증, 회전 강체 루프 |
| `HasCollision` 조기 종료 동치 | 없음 | 86,794건 불일치 0, 0.6.0 대비 시간 |
| 좌표 크기 의존성 | 라운드 2: double 경로의 거리만 | cast가 좌표 크기에 비례해 나빠짐(P1) |
| net10.0 | 미검증 (로컬 SDK 9) | SDK 10으로 실행 |
| 0.6.0 -> 0.7.0 출력 비교 | - | 1401줄 diff 0 |
