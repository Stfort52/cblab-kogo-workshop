# KOGO Workshop: Single-cell RNA-seq Analysis with R

## Introduction

This is designed to serve about ~60 users for a KOGO RNA-seq workshop.
Check [this link](https://www.kogo-edu.or.kr/workshop/info/7)

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

If you really don't have one, comment out the `caddy` service and uncomment the following to the `jupyterhub` service.

```yaml
services:
  jupyterhub:
    ports:
      - "8000:8000"
```

This will fall back to HTTP, so I recommend that you don't never handle sensitive data there!

### 4. Start the JupyterHub server

Simply run:

```bash
docker compose up --build
```
