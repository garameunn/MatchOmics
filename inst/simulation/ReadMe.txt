## ============================================================
##  MatchOmics Simulation – ReadMe
##  Last updated: 2026-04
## ============================================================

## 서버 환경
##   경로   : /home2/nekim/Matching2026/
##   R 버전 : R 4.4  (/usr/local/R-4.4/bin/Rscript)
##            ※ 'Rscript' 명령은 R 3.x 이므로 반드시 풀 경로 사용
##   구조   :
##     code/    ← R 스크립트 전체
##     data/    ← otulist, indiclist, psdifflist, seed, batchinfo
##     results/ ← T1/, Power/, Organized/, Diagnostics/

## 파일별 설명
##   00_library_setup.R          CRAN + Bioconductor 패키지 설치/로드
##   01_PreliminaryRobject.R     서버에서 실행 불필요 (data/ 이미 존재)
##   02_simulation_generator.R   핵심 시뮬레이션 함수 (수정본)
##   99_utils.R                  공용 유틸 함수 (수정본)
##   03_Type1Error_Calculating_server.R   T1 error 시뮬레이션
##   04_Power_Calculating_server.R        Power 시뮬레이션
##   05_Bias_EffectSize_Calculating_server.R  bias/effectsize
##   06_Final_Stat_Organizing_server.R    결과 정리 → xlsx
##   07_Diagnostics_server.R             SMD, Neff, Round2 진단

## 주요 변경 이력
##   - 매칭 알고리즘: with replacement → 1라운드 without replacement
##                    + 2라운드 구제 매칭 (two_round_match)
##   - 가중치 버그 수정: mm = m/sum(m) → m/mean(m)  (GEE sum(w)=N 보장)
##   - 평가 지표 추가: SE 별도 저장, Neff, SMD/Love plot, Round2 진단
##   - S1: nonlin=inv 단일 (comp 제거)
##   - 외부 loop 코어: --cores 3
##   - 내부 per_taxa loop: mc.cores=3 (코드 내 고정)
##   - 총 실효 코어: 03 × 9, 04 × 9 (한 노드 동시 실행 시 ~18 코어)

## 시뮬레이션 설계
##   시나리오  : S1 / S2 / S3
##   샘플 크기 : 200 / 500 / 1000
##   유병률    : 0.1 / 0.2 / 0.3
##   캘리퍼    : NULL / 0.3 / 0.2 / 0.1 / 0.05  (5종)
##   상관구조  : independence / exchangeable  (GEE)
##   반복 수   : T1 error = 5000회, Power = 1000회

## 실행 방법
##
##  방식 1 (한 노드):
##    bash run_simulation.sh all
##      → screen t1all: 03 전체 (3코어)
##      → screen power: 04 전체 (3코어)
##
##  방식 2 (노드 분리):
##    bash run_simulation.sh node1   → 03 S1,S2
##    bash run_simulation.sh node2   → 03 S3 + 04 전체
##
##  공통 (03+04 완료 후):
##    bash run_simulation.sh post    → 05 → 06 → 07 순차
##
##  개별 실행:
##    /usr/local/R-4.4/bin/Rscript 03_Type1Error_Calculating_server.R --cores 3
##    /usr/local/R-4.4/bin/Rscript 03_Type1Error_Calculating_server.R --cores 3 --scenario S1,S2
##    /usr/local/R-4.4/bin/Rscript 04_Power_Calculating_server.R --cores 3
##
##  모니터링:
##    screen -r t1all
##    tail -f /home2/nekim/Matching2026/results/log_03.txt
##    tail -f /home2/nekim/Matching2026/results/log_04.txt

## 결과 파일 구조
##   results/T1/      NSAMP{n}_REP5000_RATE{r}_inv_{S}.RData
##                    NSAMP{n}_REP5000_RATE{r}_inv_{S}_SMD.RData
##   results/Power/   NSAMP{n}_REP1000_RATE{r}_inv_{S}.RData
##   results/Organized/  T1_organized.RData, Power_organized.RData
##                        T1_combined.xlsx, Power_combined.xlsx
##                        Beta_SE_combined.xlsx, RelBias_combined.xlsx
##                        Nrow_Neff_combined.xlsx
##   results/Diagnostics/ SMD_combined.xlsx, Neff_summary.xlsx
##                         Round2_comparison.xlsx
##                         LovePlot_*.png, Neff_barplot.png
