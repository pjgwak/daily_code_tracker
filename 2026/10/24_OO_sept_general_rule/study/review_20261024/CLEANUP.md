# 폴더 정리 기록

현재 문서와 그림의 검토 위치를 10월 `review/`로 통합했다. 기존 빈별 Markdown 문서와 목차를 이동한 뒤 선택된 전체 실행의 수치·그림·설정으로 갱신했다. 23일 폴더의 과거 결과 설명도 최신 결과로 교체했다.

| 항목 | 정리 전 | 정리 후 |
|---|---:|---:|
| `23_OO_simple_rules_2p1/results/`의 실행·bundle 항목 | 335 | 4 |
| `23_OO_simple_rules_2p1/study/`의 항목 | 410 | 2 |
| 이력 보관으로 이동한 상위 항목 | — | 754 |
| 이동한 파일 | — | 81,289 |
| 이동한 파일 크기 | — | 약 2.55 GiB |
| 삭제한 과학적 결과 파일 | — | 0 |

현재 결과·NP prefit·DATA prefit·PDF 재생 검증의 네 결과 폴더, 입력 ROOT, 기본 prefit ROOT, 소스와 헤더, calibration window map은 23일 폴더에 남겼다. 선택한 campaign 내부에는 20개 초기값 ROOT와 범위 사본이 있어 과거 스터디 폴더 없이 재현할 수 있다.

중간 실행·그림·스크립트는 [이력 보관 폴더](../archive/23_OO_intermediate_20261024/README.md)로 이동했다. 개별 이동은 같은 파일시스템의 rename으로 수행하고 inode를 대조했다. 원본 campaign manifest와 ROOT 내용은 수정하지 않았다. 과거 manifest에 기록된 원래 경로는 [이동 계획](audit/cleanup_plan.json)과 [실제 이동 기록](audit/cleanup_completed.json)으로 추적한다.

리뷰 자산의 이동 출처는 [relocation.json](audit/relocation.json)에 기록했다. 그림·문서의 현재 해시는 `audit/review_verification.json`에 기록한다. 과거 manifest에 포함된 당시 경로와 현재 문서의 활성 링크는 구별한다.

## 현재 작업 절차

[23일 폴더](../23_OO_simple_rules_2p1/README.md)에서 `workflow/verify_review.py`로 입력·ROOT·그림·링크를 확인하고, `workflow/refresh_review.py`로 빈별 문서를 갱신한다. `workflow/reproduce_fit.py`는 선택된 실행 내부의 코드·규칙·초기값·범위 사본을 사용해 별도 campaign으로 재현한다.

보관된 중간 스터디 스크립트는 당시 디렉터리 구조를 전제로 하므로 현재 workflow와 구분한다. 현재 결과를 다른 실행으로 자동 교체하는 동작은 추가하지 않았다.

## 정리 후 검사

현재 입력 126개·override prefit 61개·최종 ROOT 20개·초기값 20개·범위 20개의 해시를 대조했다. 빈별 문서 20개와 활성 Markdown 링크를 확인했다. 정리된 구조에서 20빈의 저장 PDF를 읽기 전용으로 다시 재생해 이전 PDF 표본과 대조했으며 최대 상대 차이는 0이었다. fit을 다시 수행하지 않았다. [검증 기록](audit/after_cleanup_pdf_replay.json).
