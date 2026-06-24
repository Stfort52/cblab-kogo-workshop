import os
c = get_config()



# 1. Authenticator
c.JupyterHub.authenticator_class = "shared-password"
c.SharedPasswordAuthenticator.user_password = ""
c.SharedPasswordAuthenticator.admin_password = ""
c.Authenticator.allowed_users = {f"edu{i:02d}" for i in range(1, 61)}
c.Authenticator.admin_users = {f"edu{i:02d}" for i in range(1, 5)}

# 2. Spawner setup
c.JupyterHub.spawner_class = "dockerspawner.DockerSpawner"
# The image to spawn for users
c.DockerSpawner.image = "kogo26-jupyterlab:latest"
c.DockerSpawner.remove = True
# This allows the containers to run usermod to change UID
c.DockerSpawner.extra_create_kwargs.update({"user": "root"})  
# Healthcheck
c.DockerSpawner.extra_create_kwargs.update({
    'healthcheck': {
        'Test': ['CMD-SHELL', 'curl -f http://localhost:8888/api || exit 1'],
        'Interval': 60000000000,   # 60s
        'Timeout': 10000000000,    # 10s
        'Retries': 3,
        'StartPeriod': 15000000000 # 15s
    }
})


c.DockerSpawner.volumes = {
    "jupyterhub-user-{username}": "/home/jovyan/work",
    "/BiO/data": {
        "bind": "/home/jovyan/data",
        "mode": "ro"
    }
}
c.Spawner.mem_limit = "10G"
c.Spawner.cpu_limit = 2.0

def uid_remap(spawner):
    """Remap the UID of the user inside the container to match the host."""
    username = spawner.user.name
    # Extract the numeric part of the username (e.g., "edu01" -> 1)
    user_number = int(username[3:])
    # Calculate the new UID based on the user number
    new_uid = 1000 + user_number  # Adjust this base UID as needed
    spawner.environment["NB_UID"] = str(new_uid)

c.DockerSpawner.pre_spawn_hook = uid_remap

# 3. Networking (The crucial part for DooD)
# Get the network name from the environment variable (set in docker-compose)
network_name = os.environ.get("DOCKER_NETWORK_NAME", "jupyterhub-network")
c.DockerSpawner.network_name = network_name

# Set the cookie age to 3 days (default is 14 days)
c.JupyterHub.cookie_max_age_days = 3

# Tell the spawned containers to connect to this network
c.DockerSpawner.extra_host_config = { "network_mode": network_name }

# Because the hub is inside a container, spawned containers need to know 
# how to reach it. We use the docker-compose service name "jupyterhub"
c.JupyterHub.hub_ip = "jupyterhub"
c.JupyterHub.hub_port = 8081

# 4. Container Management
# Automatically remove stopped user containers
c.JupyterHub.cleanup_servers = True
c.JupyterHub.concurrent_spawn_limit = 40

# 5. Miscellaneous
c.JupyterHub.logo_file = "/srv/jupyterhub/logo.png"