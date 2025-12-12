# 목적
- 논문 `High-Q photonic nanocavity in a two-dimensional photonic crystal`(Akahane et al., Nature 2003) 구조를 MEEP 기반 FDTD로 재현하고, L3 결함의 gentle confinement 효과를 수치로 확인한다.

# 주요 참고 파라미터 (논문)
- 삼각 격자: 격자상수 `a = 0.42 μm`, 슬랩 두께 `T = 0.6a (~0.25 μm)`, 홀 반지름 `R = 0.29a (~0.12 μm)`.
- L3 결함: 3개 홀 제거 후 양끝 홀을 외측 변위(예: 첫 이웃 ±0.15a, 두 번째 ±0.05a), 필요 시 에지 홀 반경 축소.
- 목표 성능: `Q ≈ 4.5×10^4` 이상, 모드부피 `V ≈ (λ/n)^3`, 넓은 FSR.

# 작업 범위 (이번 단계)
- MEEP C++ API로 3D FDTD 셋업(PML, 소스, 모니터, Harminv) 구성.
- 기하 생성: 삼각 격자 파라미터화 + L3 결함 + 에지 변위/반경 조정 함수 구현.
- 밴드/주파수 설정: MPB 또는 간단한 스캔으로 TE 밴드갭 위치 확인 후 공진 탐색 대역 결정.
- Q 추출: 링다운 신호를 Harminv로 피팅, Poynting 박스에서 방사 손실 채널 분리.
- 후처리: 필드 스냅샷, k-공간 분석(필요 시 FFT), 파라미터 스윕 결과 정리.

# 현재 상태/환경
- MEEP 1.32.0-beta를 `~/.local`에 빌드·설치 완료, `meep.pc` 확인됨.
- `PKG_CONFIG_PATH`, `LD_LIBRARY_PATH`에 `~/.local` 추가됨, CMake RPATH도 `~/.local/lib`로 설정.
- C++ 바이너리에서 `#include <meep.hpp>` 및 `meep::linspace` 호출로 링크 검증 완료.

# 구체적 TODO
- [ ] MEEP C++ 시뮬레이션 스켈레톤: PML 두께, 시간 스텝, 해상도(픽셀/a) 기본값 정의.
- [ ] 기하 생성기: a, r, t, Δx1/Δx2, r_edge를 입력으로 L3 구조 빌드.
- [ ] 소스/모니터 배치: 가우시안 펄스(Ex/Ey) + Harminv 모니터 + 필드 스냅샷(중앙 평면).
- [ ] Q 계산 파이프: Harminv 결과로 Q 계산, Poynting 플럭스 박스에서 z-방사/측방 분리.
- [ ] 파라미터 스윕: Δx1/Δx2, r_edge를 소규모 그리드로 스윕하는 래퍼 작성.
- [ ] 시각화: Python 또는 C++에서 필드/플럭스/스펙트럼 저장 및 플롯 스크립트 준비.

# 모듈화 제안 (클래스 단위)
- `LatticeGeometry`: 삼각 격자/L3 결함/홀 변위·반경 설정 → MEEP geometry list 생성.
- `SimulationConfig`: 격자상수/해상도/PML/시간스텝/블록 크기 등 공통 파라미터 관리, `meep::structure`/`meep::grid_volume` 초기화.
- `SourceManager`: 가우시안 펄스/편광 설정, 소스 배치 및 초기화.
- `MonitorSuite`: Harminv 링다운, Poynting 플럭스 박스, 필드 스냅샷 설정 + 데이터 수집.
- `QAnalyzer`: Harminv 결과로 공진 주파수·Q 계산, 플럭스 분해(z/측방), 모드부피 계산.
- `SweepRunner`(선택): Δx1/Δx2/r_edge 파라미터 스윕 실행 및 결과 저장(CSV/JSON).
# 제외/보류
- MFEM 기반 FEA 해석은 보류(필요 시 후속 단계에서 곡면/복잡 형상 정밀 해석용으로 별도 분기).
