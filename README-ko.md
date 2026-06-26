# KOGO 워크숍: R을 이용한 단일 세포 RNA-seq 데이터 분석

[EN](./README.md) | KO

Single-cell RNA-seq Analysis with R

## 강좌 소개

이 프로젝트는 [한국유전체학회](https://www.kogo.or.kr/) 통계유전학워크샵의
`[단일세포전사체 기초 #1] 개념과 분석의 이해` 강좌를 위하여 만들어졌습니다.
2026년도 내용은 [이 링크](https://www.kogo-edu.or.kr/workshop/lectureInfo/7/425)에서 확인하실 수 있습니다.

### 강좌 구조

본 강좌에서는 DooD(Docker-out-of-Docker) 방식으로 JupyterLab 클라이언트를 제공하는 [JupyterHub](https://jupyter.org/hub) 서버를 기반으로 구성되어 있습니다.

## 설치하기

### 1. 저장소 가져오기

```bash
git clone https://github.com/Stfort52/cblab-kogo-workshop.git
```

### 2. 수강생 도커 이미지 (jupyterlab) 빌드하기

```bash
cd client
docker build -t kogo-workshop-client .
```

별도로 사용하고자 하는 이미지가 있으시다면 뛰어넘으셔도 무방합니다.

### 3. 강좌 설정 수정하기

#### 기본적인 설정

프로젝트의 `.env` 파일을 열어보시면 설정 항목을 찾으실 수 있습니다.
이 설정들이 런타임에 환경 변수로서 서버에 주입되어 반영되게 됩니다.

```bash
DATA_PATH=/BiO/data
IMAGE_NAME=kogo-workshop-client
N_USERS=60
N_ADMIN=4
MEM_LIMIT=8G
CPU_LIMIT=2.0
SERVER_ADDR=your-server-domain
```

- `DATA_PATH`: 워크샵 진행에 필요한 데이터 파일들의 *호스트* 경로입니다. 이 경로들이 각 수강생의 `$HOME/data` 디렉토리에 읽기 전용으로 마운트됩니다.
- `IMAGE_NAME`: 수강생들이 사용하게 될 도커 이미지입니다. 꼭 주피터랩 기반의 이미지를 사용하시길 바랍니다. 공식 주피터랩 이미지를 베이스로 만드셔도 좋고 ([기본](https://quay.io/repository/jupyter/minimal-notebook) 혹은 [R 기반으로](https://quay.io/repository/jupyter/r-notebook)) 직접 만드셔도 좋습니다.
- `N_USERS`: 관리자(조교)와 일반 사용자(수강생)을 포함한 전체 사용자 숫자입니다. 예를 들어 `N_USERS=60`은 `edu01`부터 `edu60`까지 60개의 계정을 만듭니다.
- `N_ADMIN`: 관리자 계정의 숫자입니다. 전체 계정 중 `N_ADMIN`개의 계정을 관리자 계정으로 만들게 됩니다. 예를 들어, `N_ADMIN=4`는 `edu01`부터 `edu04` 까지 관리자 권한을 부여합니다.
- `MEM_LIMIT`: 각각의 수강생들의 컨테이너가 점유할 수 있는 메모리의 상한입니다. [여기](https://docs.docker.com/engine/containers/resource_constraints/#limit-a-containers-access-to-memory)를 참조하여 자세히 알아볼 수 있습니다.
- `CPU_LIMIT`: 각각의 수강생들의 컨테이너가 동시에 점유할 수 있는 CPU 연산량의 상한입니다. 소숫점을 사용할 수 있는데, 이를테면 2.4는 최대 240%의 CPU 사용량을 허용합니다. 이때 CPU 3개가 80%로 돌아가는지 4개가 60%로 돌아가는지는 신경쓰지 않습니다. 자세히는 [여기](https://docs.docker.com/engine/containers/resource_constraints/#cpu)를 읽어보세요. Docker(와 [CFS](https://docs.kernel.org/scheduler/sched-design-CFS.html))는 기본적으로 CPU 사용량을 *논리적* 코어로 측정한다는 것을 명심하시길 바랍니다 .
- `SERVER_ADDR`: 서버를 올릴 도메인입니다. 이 주소로 [`caddy`](https://hub.docker.com/_/caddy)가 HTTPS 설정을 진행할 것입니다.

#### 인증 설정

현장에서 100% 대면으로 진행되는 1일 정도로 짧은 워크샵의 경우에는 공유 비밀번호를 사용하는 편이 편리합니다. 비밀번호를 `pass.json` 파일에서 지정하면, 런타임에 서버에 주입되어 설정됩니다.

```json
{
    "ADMIN_PASS": "AreYouTalkingAboutTheIncidentWhereChaeSooBinAndKooJaWookWereWalkingHandInHand",
    "USER_PASS": "Matgukno"
}
```

참고로, `JupyterHub`는 기본적으로 일반 사용자에 8자, 관리자에 24자의 비밀번호 글자수 하한을 적용하고 있습니다.

#### 저는 도메인이 없어요

정말 도메인이 없으면 `docker-compose.yaml`에서 `caddy` 서비스를 제거하고 `jupyterhub` 서비스에서 주석 처리된 `ports` 설정을 사용하실 수 있습니다.

```yaml
services:
  jupyterhub:
    ports:
      - "8000:8000"
```

그리고 도메인 대신 `<서버 IP 주소>:8000`을 사용하시면 됩니다. 다만 이러면 HTTP를 사용하게 되니 민감한 정보는 다루지 않는 편이 좋을지도 모릅니다.

### 4. 서버 시작하기

간단하게:

```bash
docker compose up --build
```
