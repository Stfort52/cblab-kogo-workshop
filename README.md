# KOGO Workshop: Single-cell RNA-seq Analysis with R

EN | [KO](./README-ko.md)

## Introduction

This is designed to serve about ~60 users for `[scRNA-seq basics #1] Understanding concepts and analysis` class
in a [KOGO](https://www.kogo.or.kr) Statistical Genetics and Genomics workshop.
Check [this link](https://www.kogo-edu.or.kr/workshop/lectureInfo/7/425) for the 2026 run.

### Structure

This is a [JupyterHub](https://jupyter.org/hub) server serving [JupyterLab](https://jupyter.org/) clients in a Docker-out-of-Docker (DooD) manner.

## Setup Instructions

### 1. Clone the repository

```bash
git clone https://github.com/CB-postech/cblab-kogo-workshop.git
```

### 2. Build the client (jupyterlab) Docker image

```bash
cd client
docker build -t kogo-workshop-client .
```

You can skip this step if you have your own custom image.

### 3. Edit the Configuration

#### General settings

You will find the configs in the `.env` file.
This contains the environment variables which gets injected in runtime, controlling the server.

```bash
DATA_PATH=/BiO/data
IMAGE_NAME=kogo-workshop-client
N_USERS=60
N_ADMIN=4
MEM_LIMIT=8G
CPU_LIMIT=2.0
THR_LIMIT=2
SERVER_ADDR=your-server-domain
```

- `DATA_PATH`: This is your *host* path with needed data files, to be read-only mounted to every client's `$HOME/data`.
- `IMAGE_NAME`: This is the docker image to be used for clients. Make sure to use a jupyterlab image, regardless if you base a jupyterlab image (like [this](https://quay.io/repository/jupyter/minimal-notebook) or [this](https://quay.io/repository/jupyter/r-notebook)) or build yourself.
- `N_USERS`: The number of total users, including normal users and admin. So 60 creates users from `edu01` to `edu60`.
- `N_ADMIN`: The number of admin users. This sets the first `N_ADMIN` users as admin, so 4 grants `edu01` to `edu04` admin privilages.
- `MEM_LIMIT`: The maximum amount of memory each client container can use. "10G" grants 10 gigabytes, more on [this link](https://docs.docker.com/engine/containers/resource_constraints/#limit-a-containers-access-to-memory).
- `CPU_LIMIT`: The peak CPU utilization for each client container. Can be fractional, so 2.4 grants at most 240% cpu utilization, regardless of 3 CPUs running at 80% or 4 cpus at 60%. More on [this link](https://docs.docker.com/engine/containers/resource_constraints/#cpu). Keep in mind that Docker (and the [CFS](https://docs.kernel.org/scheduler/sched-design-CFS.html) behind it) measures CPUs in *logical* cores.
- `THR_LIMIT`: The number of *threads* each R kernel will use. With this unset, some libraries might distribute the 200% `CPU_LIMIT` across 100 threads, resulting in suboptimal performance. Unlike `MEM_LIMIT` and `CPU_LIMIT`, this is not enforced but suggested to R via environment variable. We suggest `THR_LIMIT` * `# of concurrent kernels` ~ `CPU_LIMIT`.
- `SERVER_ADDR`: The domain for your server, used by [`caddy`](https://hub.docker.com/_/caddy) for HTTPS setup.

#### Authentication

For a short-lived, on-site workshop, a shared password authenticator is usually sufficient. You can set the shared password in the `pass.json` file. This will get injected in runtime too.

```json
{
    "ADMIN_PASS": "superduperlongpasswordthatyoucannotremember",
    "USER_PASS": "arelativelyshorterpassword"
}
```

Keep in mind that `JupyterHub` employs a lower bound of 24 and 8 characters for admin/user password length, respectively.

#### I don't have a domain

If you really don't have one, comment out the `caddy` service and uncomment the following to the `jupyterhub` service from the `docker-compose.yaml`.

```yaml
services:
  jupyterhub:
    ports:
      - "8000:8000"
```

This will fall back to HTTP, so I recommend that you never handle sensitive data there!

#### (Optional) Add a logo

You can add a logo to the top left of the JupyterHub page by changing the placeholder at `server/logo.png`. Why not use the [KOGO logo](https://www.kogo.or.kr/html/user/core/view/reaction/header/357/inc/images/logo.png)?

### 4. Start the JupyterHub server

Simply run:

```bash
docker compose up --build
```
