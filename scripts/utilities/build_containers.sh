#!/bin/bash

# Build Docker containers for omniCADD pipeline
# Usage: ./build_containers.sh [container_name] or ./build_containers.sh all

set -e

DOCKER_DIR="workflow/docker"
REGISTRY_PREFIX="juliahoglund"  # Change this to your Docker Hub username

build_container() {
    local container_name=$1
    local dockerfile_path="${DOCKER_DIR}/${container_name}/Dockerfile"
    
    if [[ ! -f "$dockerfile_path" ]]; then
        echo "Error: Dockerfile not found at $dockerfile_path"
        return 1
    fi
    
    echo "Building ${container_name} container..."
    cd "${DOCKER_DIR}/${container_name}"
    
    # Special handling for renv container (needs config context)
    if [[ "$container_name" == "renv" ]]; then
        cd ../../..
        docker build -f workflow/docker/renv/Dockerfile -t "${REGISTRY_PREFIX}/${container_name}:latest" --build-arg BUILDKIT_INLINE_CACHE=1 .
    else
        docker build -t "${REGISTRY_PREFIX}/${container_name}:latest" --build-arg BUILDKIT_INLINE_CACHE=1 .
        cd ../../..
    fi
    
    echo "✅ Successfully built ${REGISTRY_PREFIX}/${container_name}:latest"
}

push_container() {
    local container_name=$1
    echo "Pushing ${container_name} to Docker Hub..."
    docker push "${REGISTRY_PREFIX}/${container_name}:latest"
    echo "✅ Successfully pushed ${REGISTRY_PREFIX}/${container_name}:latest"
}

# Available containers
CONTAINERS=("mafTools" "sift4g" "gerp" "phast" "renv" "augustus" "snpeff" "vep")

case "${1:-all}" in
    "all")
        echo "Building all containers..."
        for container in "${CONTAINERS[@]}"; do
            build_container "$container"
        done
        echo ""
        echo "All containers built successfully!"
        echo "To push to Docker Hub, run: ./build_containers.sh push"
        ;;
    "push")
        echo "Pushing all containers to Docker Hub..."
        for container in "${CONTAINERS[@]}"; do
            push_container "$container"
        done
        echo "All containers pushed successfully!"
        ;;
    "mafTools"|"sift4g"|"gerp"|"phast"|"renv"|"augustus"|"snpeff"|"vep")
        build_container "$1"
        ;;
    *)
        echo "Usage: $0 [all|push|container_name]"
        echo ""
        echo "Available commands:"
        echo "  all      - Build all containers (default)"
        echo "  push     - Push all containers to Docker Hub"
        echo ""
        echo "Available containers:"
        echo "  mafTools - Build only mafTools container"
        echo "  sift4g   - Build only SIFT4G container"
        echo "  gerp     - Build only GERP container"
        echo "  phast    - Build only PHAST container"
        echo "  renv     - Build only R environment container"
        echo "  augustus - Build only Augustus container"
        echo "  snpeff   - Build only SnpEff container"
        echo "  vep      - Build only VEP container"
        exit 1
        ;;
esac