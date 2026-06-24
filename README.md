# KOGO Workshop: Single-cell RNA-seq Analysis with R

## Introduction

This is designed to serve about ~60 users for a KOGO RNA-seq workshop.
Check [this link](https://www.kogo-edu.or.kr/workshop/info/7)

### Structure

This is a [JupyterHub](https://jupyter.org/hub) server serving [JupyterLab](https://jupyter.org/) clients in a Docker-out-of-Docker (DooD) manner.

## Setup Instructions

### 1. Clone the repository

```bash
git clone (url will go here later)
```

### 2. Build the client (jupyterlab) Docker image

```bash
cd client
docker build -t kogo-workshop-client .
```

### 3. Edit the JupyterHub configuration file

You will find the `jupyterhub_config.py` file in the `server` directory.
You might find the following parts useful to edit:

#### 1. Authenticator

For a short-lived, on-site workshop, a shared password authenticator is usually sufficient. You can set the shared password in the `jupyterhub_config.py` file:

```python
c.SharedPasswordAuthenticator.user_password = "your_password_here" # more than 8 characters
c.SharedPasswordAuthenticator.admin_password = "your_admin_password_here" # more than 24 characters
```

Also, change the number of users and admins if needed:

```python
c.Authenticator.allowed_users = {f"edu{i:02d}" for i in range(1, 61)}
c.Authenticator.admin_users = {f"edu{i:02d}" for i in range(1, 5)}
```

If you want to use more advanced authentication methods, you can refer to the [JupyterHub documentation](https://jupyterhub.readthedocs.io/en/stable/getting-started/authenticators-users-basics.html). A PAM-based authenticator is a good option for a more secure setup, allowing users to log in with their own credentials.

#### 2. Spawner setup

If you want to use a different Docker image for the workshop, you can change the `image` parameter in the `DockerSpawner` configuration. Just be sure to build the Docker image before starting the JupyterHub server.

```python
c.DockerSpawner.image = "your-image-name:tag"
```

Also, you can set resource limits for the spawned containers. For example, to limit the CPU and memory usage:

```python
c.DockerSpawner.cpu_limit = 2.0  # Limit to 2 CPU cores
c.DockerSpawner.mem_limit = "4G"  # Limit to 4 GB of memory
```

#### 3. Volume directory

The data directory is set to `/BiO/data` by default.
Change this bind mount to point appropirate directory.

```py
c.DockerSpawner.volumes = {
    "jupyterhub-user-{username}": "/home/jovyan/work",
    "/BiO/data": {
        "bind": "/home/jovyan/data",
        "mode": "ro"
    }
}
```

### 4. Start the JupyterHub server

#### Configure HTTPS

Edit the caddy's command so that it points to your domain

```yaml
services:
  caddy:
    # Insert a proper domain and point it to your server
    command: caddy reverse-proxy --from <your domain here> --to  jupyterhub:8000
```

If you really don't have one, comment out the `caddy` service and uncomment the following to the `jupyterhub` service.

```yaml
services:
  jupyterhub:
    ports:
      - "8000:8000"
```

This will fall back to HTTP.

#### Launch the server

Simply run:

```bash
docker compose up --build
```
