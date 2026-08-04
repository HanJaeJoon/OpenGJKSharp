# OpenGJKSharp 0.5.0 dogfooding findings

NuGet 패키지 OpenGJKSharp 0.5.0을 net8.0 콘솔 프로젝트(`Dogfood/`)에서 실사용 검증한 기록.
1.0.0 체크리스트: (1) ISupportMappable 실용성, (2) GkPolytope/GkSimplex public 필요성, (3) CsFunction 제거 영향.

## 시나리오 완료 체크리스트

- [x] 1. 캐릭터 컨트롤러 (Capsule vs 바닥/벽 depenetration)
- [x] 2. 구 무더기 정합성
- [x] 3. 커스텀 shape 구현 체험 (RoundedBox / LineSegment / Point)
- [x] 4. API 인체공학 평가 (transform 래퍼 부재 등)
- [x] 5. 성능/할당 측정
- [x] 6. 극단 케이스
- [x] 7. public 표면 감사 (GkPolytope/GkSimplex)
- [x] 8. CsFunction 감사

---

# 최종 요약: 1.0.0 진행 가능 여부와 권고 API 변경 목록

**작성일: 2026-08-01 (시나리오 8/8 완료 후)**

## 결론: 1.0.0 진행 가능. 단, 아래 P0 2건을 1.0에 함께 넣을 것을 권고.

라이브러리의 **수치 품질은 1.0 수준으로 충분하다.** 8개 시나리오 전체에서 NaN/무한/음수 거리 0건,
대칭성 오차 0, 해석값 대비 오차 대부분 float 정밀도 한계(1e-6 이하). 1e6 대좌표, 1e-3 tiny 도형,
1e-5 미세 겹침, 정확히 접촉, 퇴화 도형(1점 hull, p0==p1 캡슐)까지 모두 안정적이었다.
심각한 버그는 1건뿐이며(깊은 곡면 관통 시 EPA 5.6% 과소평가) 물리 엔진 정상 동작 범위에서는 오차 0이다.

## 3가지 검증 목표에 대한 답

**(1) ISupportMappable이 실제 물리 워크로드에서 쓸 만한가: 예, 시그니처 변경 없이 1.0 확정 가능.**
메서드 2개(GetSupportPoint/GetCenter)로 캐릭터 컨트롤러, 커스텀 도형 3종(RoundedBox/LineSegment/Point),
pose 래퍼까지 전부 구현했고 정확도 손실이 없었다. 시그니처를 바꿀 이유는 발견되지 않았다.
다만 **direction 계약(영벡터 가능성/정규화 여부)과 GetCenter의 영향 범위는 문서화가 필요하다** -
이것 없이는 구현자가 방어적으로 짤 수밖에 없다.

**(2) GkPolytope/GkSimplex는 소비자에게 필요한가: 기능적으로 불필요 -> internal화 찬성.
단, 저할당 대체 경로를 같은 릴리스에 넣어야 한다.**
시나리오 1~6에서 단 한 번도 필요하지 않았다. 반면 직접 써 보니 **유일한 저할당 경로**였다
(cube-cube 쿼리: GkPolytope 재사용 368 B vs vertex-array 1832 B vs support-mapped 4088 B).
이걸 그냥 없애면 할당에 민감한 소비자는 0.5.0에 머물 이유가 생긴다. 게다가 이 타입들은
`Coord`가 `double[][]`인데 레이아웃([n][3])이 미문서화이고 `NVrtx`까지 settable이어서
소비자가 내부 불변식을 깰 수 있다 - public에 있으면 안 되는 형태다.

**(3) CsFunction 없이 모든 시나리오가 가능한가: 예, 단 하나 예외 - double 정밀도 입력.**
수치 동등성 완전 일치(랜덤 300쌍 maxDiff = 0), 기능은 진부분집합, 마이그레이션은 5줄 변환 헬퍼로 끝.
그러나 **CsFunction은 유일한 double 입력 경로**다. 1e6 좌표에서 1e-3 간격을 CsFunction은 해석하고
Vector3[] 경로는 0으로 뭉갠다. 제거 자체는 찬성하되 double 진입점의 존폐를 의식적으로 결정해야 한다.

## 권고 API 변경 목록 (우선순위순)

### P0 - 1.0에 반드시 함께

1. **재사용 workspace API 추가** (체크리스트 2 선결 조건)
   `OpenGjk.ComputeMinimumDistance(a, b, GjkWorkspace ws)` 형태로 내부 버퍼를 노출 없이 재사용.
   근거: 현재 support-mapped 쿼리 1회에 4088 B, EPA는 151840 B. 60Hz/1000쌍이면 초당 약 245MB 할당.
   이게 있으면 GkPolytope/GkSimplex internal화가 순수한 개선이 된다. (시나리오 5, 7)
2. **GkPolytope/GkSimplex internal화** (위 1번과 동시에)
   더불어 CsFunction의 `[Obsolete]` 메시지가 GkPolytope 오버로드를 대체 API로 안내하고 있으니
   이 메시지도 함께 수정. (시나리오 7, 8)

### P1 - 1.0에 넣는 것이 바람직

3. **CsFunction 제거 시 double 입력 경로 결정**: 현대적 double 오버로드
   (`ReadOnlySpan<double>` 기반)를 추가하거나, "공개 계약은 float 정밀도"를 README에 명시. (시나리오 8)
4. **TransformedShape(pose 래퍼) 내장**: local shape + rotation + translation.
   ConvexHull을 움직일 때 배열 재구축 대비 2.8배 빠르고 할당 2.5배 적음(실측).
   현재 pose 지원이 Box만 Quaternion을 받아 도형마다 비일관. (시나리오 1, 4)
5. **문서 보강** (전부 시나리오 1~3, 6에서 실제로 막혔던 지점):
   - `GetSupportPoint(direction)`: 영벡터 가능 여부, 정규화 보장 여부
   - `GetCenter()`: 정확성에 영향 없고 수렴 속도에만 영향인지
   - 충돌 시 `closestA`/`closestB`의 의미 (실제로는 교집합 내부의 동일한 점, 표면 위 아님)
   - 충돌 시 `ComputeMinimumDistance`가 정확히 0이 아니라 ~1e-16일 수 있음
   - README EPA 절에 contactNormal 방향(a -> b) 명시 (현재 XML 문서에만 있음)
   - 완전 동심인 경우 normal 방향은 임의

### P2 - 이후 릴리스로 미뤄도 무방

6. **깊은 곡면 관통에서 EPA 깊이 5.6% 과소평가 수정** (유일한 [버그])
   구-구 겹침 1.5까지 오차 0.00%, 1.9에서 2.42%, 동심(2.0)에서 5.62%. 박스-박스는 1.9에서도 0.00%이므로
   곡면 support 도형에서 EPA polytope 세분이 덜 수렴한 채 종료되는 것으로 보임.
   재현: `ComputeCollisionInformation(new Sphere(zero,1), new Sphere(zero,1), out _)` -> -1.88755 (기대 -2.0).
   얕은 접촉에서는 오차가 0이라 통상 물리 시뮬레이션에는 영향 없음. (시나리오 6)
7. **`HasCollision`이 `ComputeMinimumDistance`보다 2배 느린 것 개선 또는 문서화**
   분리 상태에서 2108 ns vs 1061 ns. 이름에서 기대되는 성능과 반대. (시나리오 5)
8. **퇴화 도형 정책 통일**: `Sphere(r=0)`/`Box(half=0)`은 예외로 거부하는데 `Capsule(p0==p1)`은 허용.
   "두께 0 내장 도형 미지원"으로 통일 명시. (시나리오 6)
9. **부호 없는 편의 오버로드** `TryGetPenetration(a, b, out depth, out normal)`:
   depenetration에서 반환값 부호를 뒤집는 실수를 막아줌. (시나리오 1)

## 검증 자산

`Dogfood/` 프로젝트(NuGet OpenGJKSharp 0.5.0 소비, net8.0)에 시나리오별 실행 가능한 코드가 남아 있다.
`dotnet run -c Release -- <1..8>`로 개별 재실행 가능. 라이브러리 소스는 열지 않고
패키지 동봉 README + XML 문서만으로 전부 구현했다.

---

## 2026-08-01 / 시나리오 1: 캐릭터 컨트롤러

### 수행한 것

- `Dogfood/Scenario1_CharacterController.cs`: Capsule 캐릭터(반지름 0.3, 세그먼트 반길이 0.4)가
  공중에서 낙하해 바닥 Box(상면 y=0)에 착지한 뒤 +x로 걸어가 벽 Box(전면 x=3)에 막히는
  600프레임(dt=1/60) 시뮬레이션. 매 프레임 `ComputeCollisionInformation(ISupportMappable, ISupportMappable, out normal)`으로
  penetration depth/normal을 얻어 depenetration(최대 4회 반복) + 법선 방향 속도 제거.
- 구현은 패키지 동봉 README.md와 XML 문서(IntelliSense)만 보고 진행. 라이브러리 소스는 열지 않음.

### 결과: PASS (버그 없음)

- 정지 상태 y=0.7000 (기대값 halfSeg+radius=0.7 정확 일치), 벽에 막힌 x=2.7000 (기대값 3-radius=2.7 정확 일치).
- 접촉 프레임 1079회 동안 NaN 없음, normal은 항상 단위 길이, 터널링 없음.
- 바닥+벽 동시 접촉(코너)도 반복 resolution으로 안정적으로 처리됨.
- 단일 프로브 검증: 바닥에 0.01 겹친 캡슐 -> d=-0.01000, normal=(0,-1,0). 문서 규약
  ("첫 번째 shape에서 두 번째 shape를 향함")과 정확히 일치.

### 발견

- [마찰] (체크리스트 1) **shape가 전부 불변 class라서 움직이는 물체는 매 프레임/매 depenetration 반복마다
  `new Capsule(...)` 할당이 필요.** 캐릭터 컨트롤러 600프레임에서 Capsule 생성이 수천 회 발생.
  위치만 바꾸는 방법(`WithCenter`/transform 파라미터/struct화)이 없음. 물리 루프 소비자 입장에서
  가장 먼저 느끼는 마찰. 정량 측정은 시나리오 5에서 수행 예정.
  (확인: `typeof(Capsule).IsValueType == false`, Sphere/Box 동일)
- [마찰] (체크리스트 1) **contactNormal의 방향 규약이 README에는 없음.** README EPA 절은
  "contactNormal is the contact normal"이라고만 하고 방향(a->b)은 XML 문서에만 있다.
  README 예제 코드 주석에 한 줄 추가 권장.
- [제안] (체크리스트 1) depenetration 유스케이스에서는 반환값 부호 해석(음수=관통)이 한 번 헷갈리면
  반대로 밀어버리기 쉬움. `bool TryGetPenetration(a, b, out depth, out normal)` 같은
  부호 없는 편의 오버로드가 있으면 실수 여지가 줄어듦 (필수는 아님).
- (체크리스트 2, 3) 이 시나리오는 GkPolytope/GkSimplex, CsFunction 없이 완전히 구현 가능했음.

## 2026-08-01 / 시나리오 2: 구 무더기 정합성

### 수행한 것

- `Dogfood/Scenario2_SpherePileConsistency.cs`: 랜덤 Sphere 400개(반지름 0.05~1.55, 10x10x10 공간)에서
  랜덤 쌍 798개를 뽑아 (a) d(a,b)==d(b,a) 대칭성, (b) |closestA-closestB|==distance,
  (c) closest point가 각 구 표면 위에 있는지, (d) 해석적 구-구 거리와 비교, (e) HasCollision과 distance 정합성 검사.
- `Dogfood/Repro_ClosestPoints.cs`: closest point 의미 확인용 최소 재현.

### 결과: PASS (수치 품질 매우 좋음)

- 대칭성 오차 0 (완전 일치), |closestA-closestB| vs distance 오차 최대 8.5e-7,
  해석적 거리 대비 오차 최대 1.6e-6. NaN/음수 거리 없음.
- HasCollision과 distance 판정 모순 없음 (겹침 11쌍 포함).

### 발견

- [마찰] (체크리스트 1) **겹치는 쌍에서 ComputeMinimumDistance가 정확히 0.0이 아니라 ~2.2e-16을 반환.**
  README는 "returns ... (0 when they collide)"라고 하므로 `distance == 0`으로 충돌 판정하는 소비자는 놓친다.
  실질적 해는 epsilon 비교지만, 문서에 "0에 매우 가까운 값일 수 있음"을 명시하거나 내부에서 0으로 클램프 권장.
  재현 (Repro_ClosestPoints.cs, 겹치는 두 구):
  ```
  A(center=(3.772,-3.275,-2.617), r=1.238), B(center=(3.382,-3.407,-2.993), r=0.760)
  ComputeMinimumDistance = 2.2245912190221953E-16   // 기대: 0.0
  ```
- [마찰] (체크리스트 1) **겹치는 쌍에서 out closestA/closestB의 의미가 미문서화.** 실제로는
  closestA == closestB인 교집합 내부의 임의 점이 나온다 (표면 위 아님, 위 재현에서 표면과 0.89 오차).
  분리 상태에서는 witness point가 정확하므로, "충돌 시 closest point는 의미 없음(또는 교집합 내부 공통점)"을
  XML 문서/README에 명시 권장.
- (체크리스트 2, 3) GkPolytope/GkSimplex, CsFunction 불필요.

## 2026-08-01 / 시나리오 3: 커스텀 shape 구현 체험

### 수행한 것

- `Dogfood/Scenario3_CustomShapes.cs`: README/XML 문서만 보고 ISupportMappable로
  RoundedBox(box support + r*normalize(dir)), LineSegment, PointShape 3종을 직접 구현.
  point-point(동일점 포함), segment-point(내부/끝점), segment-segment(skew),
  roundedbox-sphere(축/코너 방향), roundedbox-roundedbox 관통(EPA), segment 관통, 정확히 접촉까지
  해석적 기대값과 비교 검증.

### 결과: PASS (12/12 검증 통과)

- 모든 거리/깊이/normal이 해석값과 일치 (오차 <= 1e-3, 대부분 float 정밀도 수준).
- 두께 0 도형(점, 선분)과 동일점 퇴화 케이스에서도 NaN/행업 없음.
- README의 Point 예제 덕에 진입 장벽이 낮았음. 인터페이스가 메서드 2개뿐이라 구현 자체는 5분 거리.

### 발견

- [마찰] (체크리스트 1) **GetSupportPoint의 direction 계약이 미문서화.** (a) 영벡터가 들어올 수 있는지,
  (b) 정규화되어 들어오는지 문서에 없음. RoundedBox처럼 normalize(direction)이 필요한 shape는
  영벡터 방어 코드를 넣어야 할지 판단 불가 (이번 실행에서 영벡터 호출 0회였지만 계약이 없으면 방어적으로 짤 수밖에 없음).
  XML 문서에 "direction은 영벡터가 아닐 수 있음/없음, 정규화 보장 안 됨" 명시 권장.
- [마찰] (체크리스트 1) **GetCenter 문서("An interior point, used as the initial guess for GJK")가
  결과에 미치는 영향을 알기 어려움.** 부정확한 center를 주면 틀린 결과가 나오는지, 수렴만 느려지는지
  소비자가 알 수 없음. "정확성에는 영향 없고 수렴 속도에만 영향" 같은 보장 명시 권장 (그게 사실이라면).
- (체크리스트 2, 3) GkPolytope/GkSimplex, CsFunction 불필요. ISupportMappable만으로 커스텀 도형 완결.

## 2026-08-01 / 시나리오 4: API 인체공학 평가

### 수행한 것

- `Dogfood/Scenario4_Ergonomics.cs`: "pose(회전+이동) + shape" 물리 스타일 소비 패턴 평가.
  (A) 라이브러리에 없는 TransformedShape 래퍼(world support(d) = R*localSupport(R^-1*d)+t)를
  직접 구현해 내장 oriented Box와 랜덤 pose 200개에서 교차 검증,
  (B) Capsule 회전(수동 endpoint 회전 vs 래퍼), (C) 움직이는 ConvexHull(64점) 2000프레임을
  "매 프레임 정점 배열 재구축" vs "로컬 hull + SetPose 래퍼"로 비교.

### 결과: PASS (수치는 전부 정합)

- 직접 만든 래퍼와 내장 oriented Box 거리 차이 최대 9.2e-7. 회전 캡슐도 해석값(2.5) 일치.
- 움직이는 ConvexHull 2000프레임: 재구축 방식 28.0ms / 8690 B/frame vs 래퍼 방식 10.0ms / 3490 B/frame.
  래퍼가 2.8배 빠르고 할당 2.5배 적음 (그런데도 프레임당 ~3.5KB가 남는데 이는 쿼리 내부 할당,
  시나리오 5에서 정밀 측정).

### 발견

- [마찰] (체크리스트 1) **pose 지원이 shape마다 비일관.** Box만 Orientation(Quaternion)을 받고
  Sphere/Capsule/ConvexHull은 회전 개념이 없음. Capsule 회전은 endpoint를 소비자가 직접 회전해야 하고,
  ConvexHull은 "world-space points"라서 움직이려면 배열 재구축이 기본 경로가 됨 (측정상 느리고 할당 큼).
- [제안] (체크리스트 1) **TransformedShape(pose 래퍼)를 라이브러리에 내장 권장.** 소비자가 ~30줄로
  만들 수는 있지만(위 검증처럼 정확도 손실 없음), R^-1 방향 변환 공식을 아는 사용자만 가능.
  mutable한 SetPose까지 있으면 물리 루프에서 재구축/재할당 없이 재사용 가능 - 시나리오 1의
  "매 프레임 new Capsule" 마찰도 이걸로 해소됨.
- [제안] (체크리스트 1) ISupportMappable 오버로드에 `HasCollision(a, b, precision)`이 있는 것은 좋음
  (vertex-array 버전과 대칭). XML 문서 품질도 균일했음. 다만 ConvexHull ctor의 "The array is copied"처럼
  Box/Capsule에도 스레드 안전성/불변성 한 줄이 있으면 좋겠음.
- (체크리스트 2, 3) GkPolytope/GkSimplex, CsFunction 불필요.

## 2026-08-01 / 시나리오 5: 성능/할당 측정

### 수행한 것

- `Dogfood/Scenario5_PerfAlloc.cs`: sphere-box 쿼리 각 10만 회, Stopwatch +
  `GC.GetAllocatedBytesForCurrentThread()`로 쿼리당 시간/할당 측정. Release 빌드, warmup 1000회,
  2회 실행해 재현성 확인 (할당 바이트는 회차 간 완전 동일).

### 결과 (net8.0 Release, 2회 실행 중 대표값)

```
HasCollision (separated)                          2108 ns/query,    4088 B/query
HasCollision (overlap)                             698 ns/query,    2520 B/query
ComputeMinimumDistance (separated)                1061 ns/query,    4088 B/query
ComputeMinimumDistance+witness (separated)        1059 ns/query,    4088 B/query
ComputeCollisionInformation (separated)           1112 ns/query,    4184 B/query
ComputeCollisionInformation (overlap, EPA)       18051 ns/query,  151840 B/query
new Sphere + new Box + ComputeMinimumDistance      959 ns/query,    4176 B/query
vertex-array ComputeMinimumDistance (cube-cube)   1112 ns/query,    1832 B/query
```

### 발견

- [버그 아님/마찰] (체크리스트 1) **쿼리당 할당이 물리 엔진 내부 루프 기준으로 과다.**
  분리 상태 단순 sphere-box 쿼리 1회에 약 4KB 할당. 60Hz에서 1000쌍이면 초당 약 245MB 할당 ->
  Gen0 GC를 계속 유발. 참고로 shape 객체 생성(new Sphere + new Box, 약 90B)은 전체 할당의 2%뿐이라
  **시나리오 1/4에서 지적한 "매 프레임 new Shape" 마찰은 실제 병목이 아니고, 쿼리 내부 배열 할당이 진짜 문제**임이 정량 확인됨.
  (권고: 재사용 가능한 작업 버퍼 - `GjkContext`/`QueryWorkspace` 같은 것을 소비자가 만들어 넘기는
  오버로드, 또는 내부 `ArrayPool<double>` 사용. 1.0에서 API를 굳히기 전에 결정해야 할 사항.)
- [마찰] (체크리스트 1) **EPA(겹침) 경로가 분리 경로보다 약 17배 느리고 할당은 약 37배(151KB/query).**
  깊이/normal이 필요한 depenetration은 접촉 쌍마다 이 비용을 낸다. 접촉 100쌍/프레임이면
  프레임당 약 15MB 할당. EPA polytope 버퍼 재사용이 가장 효과 큰 최적화 지점.
- [마찰] (체크리스트 1) **`HasCollision`(분리 상태)이 `ComputeMinimumDistance`보다 2배 느림
  (2108 vs 1061 ns).** 이름만 보면 boolean 조기 종료로 더 빠를 것으로 기대하는데 반대다.
  같은 할당(4088B)이므로 내부적으로 거리 계산 후 비교하는 구조로 보이는데, 그렇다면 왜 2배 느린지
  이유가 불명확(측정은 2회 재현). 소비자 기대와 어긋나므로 문서에 "HasCollision은 최적화 경로가 아님"
  명시 또는 실제 조기 종료 최적화 권장.
- [제안] (체크리스트 1) support-mapped 경로(4088B)가 vertex-array 경로(1832B)보다 할당이 2.2배 큼.
  ISupportMappable 어댑터 계층에 추가 버퍼가 있는 듯. support API를 권장 경로로 밀 것이라면
  이 격차를 줄이는 게 좋음.
- (체크리스트 2, 3) GkPolytope/GkSimplex, CsFunction 불필요. 다만 성능 튜닝 목적으로
  내부 버퍼에 접근하고 싶은 소비자에게는 GkPolytope가 유혹거리가 될 수 있음 - 시나리오 7에서 검증.

## 2026-08-01 / 시나리오 6: 극단 케이스

### 수행한 것

- `Dogfood/Scenario6_EdgeCases.cs`: (A) 동심 도형(구/구, 동일 구, 동심 박스, 박스 안 구),
  (B) 정확히 접촉(거리 0: 구-구, 박스 면-면), (C) 1e-2~1e-5 미세 겹침,
  (D) 좌표 1e6 스케일(작은 도형/큰 도형/큰 도형 관통), (E) 반지름 1e-3 tiny 도형 및 1e6에 놓인 tiny 도형,
  (F) 겹침 비율별 EPA 깊이 정확도 스윕(구-구 vs 박스-박스), (G) 퇴화 도형(반지름 0, 두께 0 박스, 1점 hull).

### 결과: NaN/무한/음수 거리 0건. 정확도 이슈 1건 발견 (아래)

- 미세 겹침: overlap 1e-5까지 깊이 오차 1.4e-8 수준 (매우 좋음).
- 정확히 접촉: 구-구/박스-박스 모두 distance = 0.000, EPA info = 0.000, normal = (1,0,0) 정상.
- 1e6 좌표: 구 gap 3 -> 3 (오차 0), 1e5 박스 gap 1e5 -> 100000 (상대오차 0), 관통 1e4 -> -10000 정확.
  대좌표에서도 안정적.
- tiny 도형(r=1e-3): gap 3e-3 -> 오차 3.2e-10. 1e6에 놓인 tiny 도형은 d=0으로 부정확하지만,
  1e6에서 float 간격이 6.25e-2라 3e-3 gap 자체가 표현 불가 - README에 명시된 입력 양자화의 결과로 예상된 동작.
- 퇴화: `Sphere(r=0)`과 `Box(halfExtents.Y=0)`은 생성 시 ArgumentException으로 거부(명확한 메시지, 좋음).
  `Capsule(p0==p1)`은 문서대로 허용되고 구로 동작(정확). 1점 `ConvexHull`도 정상 동작(d=3 정확).

### 발견

- [버그] (체크리스트 1) **깊게 겹친 곡면(구) 쌍에서 EPA 관통 깊이가 최대 5.6% 과소평가된다.**
  구-구(각 r=1) 겹침 비율 스윕 결과, 겹침 1.5까지는 오차 0.00%인데 그 이상에서 급격히 벌어진다:
  ```
  overlap   depth      relErr
    1.500    1.49997    0.00%
    1.900    1.85405    2.42%   <- 기대 1.900
    2.000    1.88755    5.62%   <- 동심(동일 구), 기대 2.000
  ```
  동심 구 r=1 안에 r=2: 깊이 2.83132 (기대 3.0, 5.6% 과소). 같은 스윕을 박스-박스로 하면
  겹침 1.9까지 오차 0.00%이므로, 다면체가 아닌 **곡면 support 도형에서 EPA polytope 세분이
  충분히 수렴하기 전에 종료**하는 것이 원인으로 보인다 (박스는 EPA가 정확히 표현 가능).
  최소 재현:
  ```csharp
  var a = new Sphere(Vector3.Zero, 1f);
  var b = new Sphere(Vector3.Zero, 1f);   // 완전 동심, 동일 반지름
  double info = OpenGjk.ComputeCollisionInformation(a, b, out _);
  // 기대: info == -2.0 (MTD 2.0)
  // 실제: info == -1.88755  (5.62% 과소평가)
  ```
  실용상 영향: 얕은 접촉(물리 엔진 정상 동작 범위)에서는 오차 0이므로 치명적이지 않다. 다만
  터널링 복구/깊은 관통 해소 시 depenetration이 부족해 여러 프레임이 걸린다. 1.0 전에
  EPA 종료 기준(tolerance/최대 반복)을 완화하거나, 곡면 도형의 정확도 한계를 문서화할 것을 권고.
- [제안] (체크리스트 1) 동심 도형처럼 EPA 초기 simplex가 원점을 포함하는 축퇴 상황에서 반환된
  normal이 임의 방향((-0.654,-0.704,0.276) 등)인데, 이는 수학적으로 불가피하다. 문서에
  "완전 동심인 경우 normal 방향은 임의" 한 줄 명시 권장.
- [제안] (체크리스트 1) `Sphere(r=0)`/`Box(half=0)` 거부는 합리적이지만 `Capsule(p0==p1)`은 허용되어
  퇴화 허용 정책이 도형마다 다르다. 1.0에서 "두께 0 도형은 지원하지 않음"으로 통일 명시 권장
  (커스텀 ISupportMappable로는 시나리오 3처럼 점/선분이 잘 동작하므로, 내장 도형만의 정책 문제).
- (체크리스트 2, 3) GkPolytope/GkSimplex, CsFunction 불필요.

## 2026-08-01 / 시나리오 7: public 표면 감사 (GkPolytope/GkSimplex)

### 수행한 것

- `Dogfood/Scenario7_PublicSurfaceAudit.cs`: 리플렉션으로 패키지 public 표면 전체 열거,
  Models 타입들의 생성/설정 가능성 확인, GkPolytope/GkSimplex 오버로드를 소비자 입장에서 실제로 호출,
  vertex-array 경로와 할당량 비교, 좌표 레이아웃 추측을 비대칭 도형으로 교차 검증.

### 결과

public 타입은 8개뿐 (Models 내부 타입 GkBody/EpaPolytope/EpaFace/EpaEdge는 이미 internal - 잘 되어 있음):

```
OpenGJKSharp.ISupportMappable       interface
OpenGJKSharp.Models.GkPolytope      class, public ctor, settable NumPoints/SIdx/Coord
OpenGJKSharp.Models.GkSimplex       class, public ctor, settable NVrtx
OpenGJKSharp.OpenGjk                static
OpenGJKSharp.Shapes.{Box,Capsule,ConvexHull,Sphere}
```

- 시나리오 1~6 전부 GkPolytope/GkSimplex 없이 구현 완료. **소비자 기능 관점에서는 불필요.**
- 그런데 직접 써 보니 **의미 있는 이점이 하나 있었다: 할당량**.
  ```
  같은 cube-cube 최소거리 쿼리, 1만 회 평균
    ComputeMinimumDistance(Vector3[], Vector3[])            1832 B/query
    ComputeMinimumDistance(GkPolytope, GkPolytope, GkSimplex)  368 B/query   <- 5배 적음
    ComputeMinimumDistance(ISupportMappable, ISupportMappable) 4088 B/query  <- 시나리오 5
  ```
  GkPolytope/GkSimplex를 재사용하면 프레임 간 버퍼가 유지되어 할당이 368B까지 떨어진다.
  현재 0.5.0에서 **저할당 경로는 이것뿐**이다.

### 발견

- [제안] (체크리스트 2) **GkPolytope/GkSimplex를 internal화하는 것은 기능적으로 타당하다** -
  시나리오 1~6에서 단 한 번도 필요하지 않았고, 8개 public 타입 중 2개가 "내부 알고리즘 상태"라는 점이
  API를 이해하기 어렵게 만든다. 특히:
  - `GkPolytope.Coord`는 `double[][]`인데 레이아웃([n][3] vs [3][n])이 문서에 없다. 이번에 [n][3]으로
    추측해 맞췄지만, 대칭 도형(cube)으로는 두 레이아웃이 구분되지 않아 비대칭 tetrahedron으로
    별도 교차 검증해야 했다. 이런 타입이 public API에 있는 것은 소비자에게 위험하다.
  - `GkSimplex.NVrtx`가 settable이고 `Coord`/`SIdx`도 settable이라, 소비자가 임의 상태를 주입해
    라이브러리 내부 불변식을 깨뜨릴 수 있다.
  - 문서 문구도 내부 구현 서술("For the first iteration this value is a guess - and this guess not
    irrelevant.")이라 소비자용이 아니다.
- [제안] (체크리스트 2) **단, internal화 전에 저할당 대체 경로를 제공해야 한다.** 위 측정대로
  GkPolytope 재사용은 vertex-array 대비 5배, support-mapped 대비 11배 할당이 적다. 이걸 그냥 없애면
  할당에 민감한 소비자(물리 엔진)는 0.5.0에 머무를 이유가 생긴다. 권고:
  1.0에서 GkPolytope/GkSimplex는 internal화하되, 시나리오 5에서 제안한 재사용 workspace
  (`OpenGjk.ComputeMinimumDistance(a, b, GjkWorkspace ws)` 형태, 내부 버퍼만 노출 없이 보유)를
  같은 릴리스에 넣을 것. 그러면 internal화가 순수한 개선이 된다.
- (체크리스트 3) CsFunction 여전히 불필요 - 시나리오 8에서 최종 감사.

## 2026-08-01 / 시나리오 8: CsFunction 감사

### 수행한 것

- `Dogfood/Scenario8_CsFunctionAudit.cs`: (A) [Obsolete] 메타데이터 확인,
  (B) 랜덤 폴리토프 쌍 300개(각 1~11점)에서 CsFunction vs ComputeMinimumDistance 수치 동등성,
  (C) 기능 비교표, (D) 입력 정밀도 차이 실측, (E) before/after 마이그레이션 샘플을 실제로 실행해 검증.

### 결과

- `[Obsolete(IsError=false)]`, 메시지에 대체 API가 명시되어 있음 (좋음). 빌드 시 CS0618 경고 정상 발생.
- **수치 동등성 완전 일치: 300쌍 maxDiff = 0.00E+000.** CsFunction은 ComputeMinimumDistance의
  기능적 진부분집합(거리만 반환, witness/depth/normal/tolerance/2D/support 도형 없음).
- 마이그레이션은 자명함: `double[3,n]` -> `Vector3[]` 변환 헬퍼 한 개(5줄)로 끝. 실행 검증에서
  before/after/헬퍼 경유 모두 3.0 동일.

### 발견

- [버그 아님/중요] (체크리스트 3) **CsFunction은 유일한 double 정밀도 입력 경로다.** 이것이
  "CsFunction으로만 가능한 것"의 유일한 사례이며, 1.0에서 제거하면 소비자가 잃는 실제 기능이다.
  최소 재현:
  ```csharp
  // 1e6 좌표에서 1e-3 간격 (float로는 표현 불가, 1e6에서 float 간격은 약 6.25e-2)
  var a = new double[3, 1]; a[0, 0] = 1e6;
  var b = new double[3, 1]; b[0, 0] = 1e6 + 1e-3;
  OpenGjk.CsFunction(1, a, 1, b);   // = 0.001000000047 (double 반올림 오차 4.75e-11, 정상)

  var fa = new[] { new Vector3(1e6f, 0, 0) };
  var fb = new[] { new Vector3((float)(1e6 + 1e-3), 0, 0) };
  OpenGjk.ComputeMinimumDistance(fa, fb);   // = 0  (간격이 완전히 소실)
  ```
  원점 근처 검증으로 double이 끝까지 보존되는 것도 확인: 거리 1/3 -> 0.3333333333333333
  (double 1/3와 비트 일치, float 1/3인 0.33333334와 구별됨). 즉 CsFunction은 이름/시그니처가
  나쁠 뿐 정밀도 면에서는 열등하지 않다.
- [제안] (체크리스트 3) 따라서 1.0에서 CsFunction 제거는 **찬성하되, double 진입점을 별도로 남길 것**을 권고.
  선택지: (a) `ComputeMinimumDistance(ReadOnlySpan<double> a, int nA, ReadOnlySpan<double> b, int nB)`
  같은 현대적 double 오버로드 추가, (b) 아무것도 안 하고 "이 라이브러리의 공개 계약은 float 정밀도"라고
  README에 명시. 대좌표(지구 규모 좌표계 등) 소비자가 있다면 (a)가 맞다.
  README에 이미 입력 양자화 언급이 있으나, "double 입력이 필요하면 X를 쓰라"는 안내는 없다.
- (체크리스트 2) CsFunction 감사 자체에는 GkPolytope/GkSimplex 불필요. 다만 [Obsolete] 메시지가
  "or ComputeMinimumDistance(GkPolytope, GkPolytope) instead"로 **1.0에서 internal화 예정인 타입을
  대체 API로 안내**하고 있다. 시나리오 7 권고대로 internal화한다면 이 메시지도 함께 고쳐야 한다.

