# OpenGJKSharp 3차 dogfooding (0.6.0-rc.1)

- 날짜: 2026-08-03
- 대상 패키지: `OpenGJKSharp 0.6.0-rc.1` (로컬 피드 `./feed`, TFM: netstandard2.0 / netstandard2.1 / net8.0)
  - net10.0은 로컬 SDK(9.0.200)가 타겟팅을 지원하지 않아 rc 패키지에서 빠졌다. 릴리스 패키지는 CI에서 4개 TFM으로 만들어진다.
- 소비자 프로젝트
  - `Dogfood3/` (net8.0) - 시나리오 1~5, 8. `dotnet run -c Release` 또는 `dotnet run -c Release -- 3`
  - `Dogfood3/Oracle.cs` - 라이브러리와 독립적인 2D 오라클 (SAT 기반 최소이동량, 선분 스윕 기반 거리)
  - `Upgrade050/` (net8.0) - 시나리오 6. `dotnet run -c Release -p:GjkVersion=0.5.0`으로 같은 소스를 0.5.0에도 빌드
  - `Upgrade050/break/` - 0.5.0의 `GkPolytope`/`GkSimplex` 코드가 0.6.0에서 어떤 진단을 받는지 확인
  - `NetFxConsumer/` (net48) - 시나리오 7
- 이번 라운드 초점: 평면/퇴화 구성 + 0.5.0 사용자 호환성

## 결론

**P0/P1 없음.** 8개 시나리오 전부 통과. 릴리스 진행에 반대할 근거를 찾지 못했다.
발견 사항은 P2 1건, friction 2건, 미검증 1건이다.

## 시나리오별 결과

### 1. 같은 평면 폴리곤 쌍 - 독립 오라클과 대조 [통과]

무작위 2D 볼록 폴리곤 400쌍(정점 3~7개)을 라이브러리와 SAT 오라클로 비교했다.

- 겹치는 289쌍의 깊이 차이 최대 `1.48e-7`, 법선 차이 최대 `1.20e-7` (float 입력 양자화 수준)
- 289쌍 전부 보고된 이동량을 적용하면 실제로 분리됨 (오라클로 재확인)
- 오라클이 "떨어져 있다"고 한 쌍에서 라이브러리가 겹침을 보고한 경우 0건
- 계획서의 재현 케이스: `Vector3[]` 같은 평면 사각형 -> `-0.250000`, 법선 `<1, 0, 0>`
- 같은 케이스의 flat `ConvexHull` 쌍 -> `-0.250000`, 법선 `<1, 0, 0>`

### 2. `Vector2[]` 2D 파이프라인 [통과]

거리/충돌/깊이/법선/witness가 한 벌로 다 되고, workspace 변형이 할당 버전과 비트 단위로 일치했다.

- 떨어진 쌍 거리 `3.000000` (오라클 일치), witness `<2, 0> <-> <5, 0>` (양쪽 표면 위)
- 겹친 쌍 `info = -0.5`, `TryGetPenetration -> depth 0.5`, 법선 `<1, 0>`
- 정확히 접촉하는 쌍: `true` + `depth 0`, 부호 비트 확인해 `-0.0` 아님

### 3. 퇴화 도형 [통과]

라운드 2까지는 이 조합에서 조용히 0이 나왔다. 이번에는 전부 기하학적으로 맞는 값이 나왔다.

| 구성 | 결과 |
| --- | --- |
| 같은 점 vs 같은 점 | depth 0 (접촉, 올바름) |
| 점이 선분 내부 | depth 1.0 = 가까운 끝까지 거리 |
| 일직선 선분 쌍 [0,2] vs [1,3] | depth 1.0 |
| 교차 선분 | depth 1.0 |
| 점이 납작한 사각형 내부 | depth 1.0, 법선 `<0, -1, 0>` |
| 선분이 사각형 모서리 위 | depth 0 (겹침이 0 측도, 올바름) |
| 선분이 사각형 내부 | depth 1.0 |
| 길이 0 `Capsule` vs 구 | 구 vs 구와 `-0.500016156`으로 소수점 9자리까지 동일 |
| 점 1개 `ConvexHull` vs `Box` | depth 1.0 |
| 선분 `ConvexHull` vs `Box` | depth 1.0 |

검증 과정에서 **내 기대값이 두 번 틀렸다**(점이 선분 내부일 때 0을 기대, 선분이 모서리 위일 때 양수를 기대).
둘 다 라이브러리가 맞았다. 퇴화 구성의 "깊이"가 직관과 다를 수 있다는 뜻이므로 아래 friction에 적었다.

### 4. 평면과 얇은 3D의 경계 [통과]

두께를 1e-3부터 1e-9까지 줄여가며 확인했다. 전 구간에서 **두께 방향 깊이**(더 작은 쪽, 즉 올바른 MTV)를
정확히 반환했고, 법선은 `<0, 0, 1>`이었다.

| 두께 | depth | 법선 |
| --- | --- | --- |
| 1e-3 | 1.000000e-3 | `<0, 0, 1>` |
| 1e-4 | 1.000000e-4 | `<0, 0, 1>` |
| 1e-5 | 1.000000e-5 | `<0, 0, 1>` |
| 1e-6 | 1.000000e-6 | `<0, 0, 1>` |
| 1e-7 | 1.000000e-7 | `<0, 0, 1>` |
| 1e-8 | 1.000000e-8 | `<0, 0, 1>` |
| 1e-9 | 1.000000e-9 | `<0, 0, 1>` |
| 0 (평면) | 0.250000 | `<1, 0, 0>` (평면 내) |

두께 0에서 값이 점프하는 것은 문서화된 의도된 동작이다(README "2D polygons and other flat geometry",
CHANGELOG Changed 항목). 평면 판정이 얇은 3D 쌍을 가로채는 회귀는 관측되지 않았다.

### 5. workspace 오용 [통과]

- 하나를 4000개 병렬 쿼리에서 공유: **3992건이 `InvalidOperationException`**, 그 외 예외 0건
- 메시지가 원인과 해법을 다 알려준다: "This GjkWorkspace is already in use by another query.
  A GjkWorkspace is not thread-safe: use one per thread instead of sharing one across threads,
  and do not start a query with the workspace of a query that is still running (for example from
  inside a custom GetSupportPoint)."
- 예외 폭풍 이후 같은 workspace로 정상 쿼리 가능 (`-0.800000`, 단일 스레드 값과 동일)
- `ThreadLocal<GjkWorkspace>`로 4000 병렬 쿼리: 값이 다른 경우 0건

라운드 2에서 "4000건 중 3618건이 조용히 틀린 값"이었던 자리다. 이제 틀린 값이 아니라 예외가 나온다.

### 6. 0.5.0 -> 0.6.0 업그레이드 [통과]

(a) **컴파일 깨짐**: `GkPolytope`/`GkSimplex`를 쓰는 0.5.0 코드는 0.5.0에서 오류 0개로 빌드되고,
0.6.0에서는 `CS0122`(보호 수준) 2건으로 실패한다. README "Migrating from 0.5.0" 절이
CS0122/CS0246과 before/after 대체 코드를 모두 담고 있어 설명되지 않는 깨짐은 없다.

(b) **값 변경**: 0.5.0 공개 API만 쓰는 소비자를 양쪽 패키지로 돌려 출력을 diff했다. 달라진 줄은 4개뿐이다.

```
-info normal                                    <1, 0, 0>     (cube 쌍)
+info normal                                    <1, -0, 0>
-info normal                                    <1, -0, 0>    (sphere/box 쌍)
+info normal                                    <1, 0, 0>
-info(concentric unit spheres)                  -1.8875492230649304
+info(concentric unit spheres)                  -1.999999959155276
-info(coplanar quads)                           0
+info(coplanar quads)                           -0.25
```

- 동심 구, 같은 평면 사각형: CHANGELOG와 README 마이그레이션 절에 있는 의도된 변경
- 법선의 `-0` / `+0` 뒤바뀜: 아래 P2 항목
- 나머지(거리, 2D 거리, witness, 박스 관통, CsFunction, 교차 캡슐)는 전부 동일

### 7. .NET Framework 소비자 복원/실행 [통과]

net48 콘솔 앱(PackageReference)에서 `System.Memory` 전이 의존성 포함 복원이 성공하고,
netstandard2.0 어셈블리가 .NET Framework 4.8.9337.0에서 정상 동작했다.

- vertex-array 거리 2, span 거리 2 (span 오버로드가 `System.Memory`를 실제로 쓰는 경로)
- `span TryGetPenetration` -> depth 0.288675
- workspace 있는 shape 쿼리 -> `-0.800000`
- 같은 평면 2D 관통 -> depth 0.250000
- 커스텀 `GetSupportPoint` 안에서 같은 workspace로 재진입 -> `InvalidOperationException` 발생 확인

### 8. 성능/할당 [통과]

40개 캡슐 바디 x ground 박스, `TransformedShape.SetPose` + `TryGetPenetration(workspace)` 물리 루프
(2000 프레임 워밍업 후 5000 프레임 측정):

- **0.01 B/frame** (라운드 2: 0~4 B/frame)
- **82~85 us/frame** (라운드 2: 109 us/frame)

## 발견 사항

### P2. 접촉 법선 성분에 `-0.0`이 나오고, 어느 성분이 그렇게 되는지가 0.5.0과 달라졌다

`ComputeCollisionInformation`이 돌려주는 법선의 성분이 `-0.0`일 수 있다. 0.5.0 대비 뒤바뀐 쌍이 있다
(cube 쌍은 `<1, 0, 0>` -> `<1, -0, 0>`, sphere/box 쌍은 `<1, -0, 0>` -> `<1, 0, 0>`).

- 수치적으로는 `-0.0 == 0.0`이라 비교/누적/`normal * depth` 결과가 달라지지 않는다
- 실제로 보이는 차이는 문자열 포맷(`"-0"`)과 성분으로 나누는 경우뿐이다
- 등급 P2로 두고 이번 릴리스에서는 고치지 않는다. 핫 경로에 정규화를 넣을 만한 이득이 없다

### friction. 퇴화 구성의 "깊이"가 직관과 어긋날 수 있다

시나리오 3을 쓰면서 두 번 잘못된 기대값을 적었다. 점이 선분 내부에 있으면 깊이는 0이 아니라
가까운 끝까지의 거리이고, 선분이 폴리곤 모서리 위에 정확히 놓이면 겹침이 0 측도라 깊이가 0이다.
둘 다 문서화된 계약(평면/직선 내 최소 이동량)에서 따라오는 값이지만, 예제 없이 읽으면 헷갈린다.
README에 퇴화 예시 한두 줄을 더 넣는 것이 1.0 전에 할 만한 개선이다.

### friction. `Parallel.For(0, n, _ => ... out _ ...)`가 컴파일되지 않는다

병렬 루프를 `_ =>`로 쓰면 `_`가 람다 매개변수(long)에 바인딩되어 `out _`이 `out long`으로 해석되고
`CS1503`이 난다. C# 언어 규칙 문제이고 라이브러리 잘못은 아니지만, 이 API 주변에서 병렬 루프를
쓰는 사람이 자주 만난다. README의 병렬 예제는 `pair =>`로 되어 있어 이 함정을 피한다.

### 미검증. packages.config 방식 복원

로컬에 `nuget.exe`가 없어 classic `packages.config` 프로젝트의 복원은 확인하지 못했다.
PackageReference 방식 net48 소비자는 확인했다(시나리오 7).

리스크 평가: 0.5.0의 netstandard2.0도 이미 `System.Numerics.Vectors` 의존성을 갖고 있었으므로,
0.6.0은 **같은 종류의 의존성이 하나 늘어나는 것**이고 새로운 문제 유형이 아니다.
README 마이그레이션 절에 PackageReference를 권장하는 문장이 들어가 있다.

## 라운드 2 미해결 항목 대조

| 항목 | 라운드 2 상태 | 라운드 3 결과 |
| --- | --- | --- |
| A. 평면 구성 관통 깊이 | 깊이 0 | 시나리오 1, 3, 4에서 오라클과 1.5e-7 이내 일치 |
| B. double 경로 `TryGetPenetration` | 없음 | 시나리오 7에서 net48까지 확인 |
| C. workspace 오용 | 3618/4000 조용히 틀림 | 3992/4000 예외, 메시지 확인 |
| D. 0.5.0 마이그레이션 안내 | CHANGELOG만 | README 절 + 시나리오 6에서 설명 누락 0건 |
| E. `dotnet pack` 함정 | 미기록 | `CLAUDE.md`에 기록 |
