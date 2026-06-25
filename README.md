# KOGO Workshop: Single-cell RNA-seq Analysis with R

## Introduction

This is designed to serve about ~60 users for a KOGO RNA-seq workshop.
Check [this link](https://www.kogo-edu.or.kr/workshop/info/7) for the 2026 run.

### Structure

This is a [JupyterHub](https://jupyter.org/hub) server serving [JupyterLab](https://jupyter.org/) clients in a Docker-out-of-Docker (DooD) manner.

## Setup Instructions

### 1. Clone the repository

```bash
git clone https://github.com/Stfort52/cblab-kogo-workshop.git
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
SERVER_ADDR=your-server-domain
```

- `DATA_PATH`: This is your *host* path with needed data files, to be read-only mounted to every client's `$HOME/data`.
- `IMAGE_NAME`: This is the docker image to be used for clients. Make sure to use a jupyterlab image, regardless if you base a jupyterlab image (like [this](https://quay.io/repository/jupyter/minimal-notebook) or [this](https://quay.io/repository/jupyter/r-notebook)) or build yourself.
- `N_USERS`: The number of total users, including normal users and admin. So 60 creates users from `edu01` to `edu60`.
- `N_ADMIN`: The number of admin users. This sets the first `N_ADMIN` users as admin, so 4 grants `edu01` to `edu04` admin privilages.
- `MEM_LIMIT`: The maximum amount of memory each client container can use. "10G" grants 10 gigabytes, more on [this link](https://docs.docker.com/engine/containers/resource_constraints/#limit-a-containers-access-to-memory).
- `CPU_LIMIT`: The peak CPU utilization for each client container. Can be fractional, so 2.4 grants at most 240% cpu utilization, regardless of 3 CPUs running at 80% or 4 cpus at 60%. More on [this link](https://docs.docker.com/engine/containers/resource_constraints/#cpu). Keep in mind that Docker (and the [CFS](https://docs.kernel.org/scheduler/sched-design-CFS.html) behind it) measures CPUs in *logical* cores.
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

### 4. Start the JupyterHub server

Simply run:

```bash
docker compose up --build
```
