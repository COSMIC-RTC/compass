#!/bin/bash
#
# Setup Python environment for COSMIC
# Supports both mamba (conda) and venv
#
# Usage:
#   ./setup_python_env.sh --mamba    # Create mamba environment
#   ./setup_python_env.sh --venv     # Create virtual environment
#

set -e  # Exit on error

# Color codes for output
RED=$'\033[0;31m'
GREEN=$'\033[0;32m'
YELLOW=$'\033[1;33m'
BLUE=$'\033[0;34m'
CYAN=$'\033[0;36m'
NC=$'\033[0m' # No Color

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Environment files
ENV_YML="$PROJECT_ROOT/compass-env.yml"
REQUIREMENTS_TXT="$PROJECT_ROOT/compass-requirements.txt"
ENV_NAME="compass"

# Function to print colored messages
print_info() {
    echo -e "${BLUE}ℹ${NC} $1"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

print_step() {
    echo -e "${CYAN}==>${NC} $1"
}

# Function to check if mamba is installed
check_mamba_installed() {
    if command -v mamba &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# Function to check if conda is installed
check_conda_installed() {
    if command -v conda &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# Function to install mamba (Miniforge)
install_mamba() {
    print_step "Installing Miniforge (mamba)..."
    
    # Detect OS and architecture
    OS=$(uname -s)
    ARCH=$(uname -m)
    
    case "$OS" in
        Linux)
            PLATFORM="Linux"
            ;;
        Darwin)
            PLATFORM="MacOSX"
            ;;
        *)
            print_error "Unsupported OS: $OS"
            exit 1
            ;;
    esac
    
    case "$ARCH" in
        x86_64)
            ARCH_SUFFIX="x86_64"
            ;;
        aarch64)
            ARCH_SUFFIX="aarch64"
            ;;
        arm64)
            ARCH_SUFFIX="arm64"
            ;;
        *)
            print_error "Unsupported architecture: $ARCH"
            exit 1
            ;;
    esac
    
    # Download and install Miniforge
    MINIFORGE_VERSION="latest"
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${PLATFORM}-${ARCH_SUFFIX}.sh"
    INSTALL_DIR="${HOME}/miniforge3"
    INSTALLER="/tmp/miniforge_installer.sh"
    
    print_info "Downloading Miniforge from $MINIFORGE_URL..."
    if ! curl -L -o "$INSTALLER" "$MINIFORGE_URL"; then
        print_error "Failed to download Miniforge"
        exit 1
    fi
    
    print_info "Installing Miniforge to $INSTALL_DIR..."
    bash "$INSTALLER" -b -p "$INSTALL_DIR"
    
    # Clean up installer
    rm -f "$INSTALLER"
    
    if [ -f "$INSTALL_DIR/bin/mamba" ]; then
        print_success "Miniforge (mamba) installed to $INSTALL_DIR"
        
        # Initialize conda for the current shell
        print_info "Initializing conda/mamba for bash..."
        "$INSTALL_DIR/bin/conda" init bash
        
        # Source the conda configuration for this session
        if [ -f "${HOME}/.bashrc" ]; then
            # shellcheck disable=SC1090
            source "${HOME}/.bashrc"
        fi
        
        print_success "Mamba installed successfully!"
        print_warning "Please restart your shell or run: source ~/.bashrc"
        print_info "Then run this script again to create the environment."
        exit 0
    else
        print_error "Failed to install Miniforge"
        exit 1
    fi
}

# Function to setup mamba environment
setup_mamba_env() {
    print_step "Setting up mamba environment..."
    
    # Check if mamba is available
    MAMBA_CMD=""
    if check_mamba_installed; then
        MAMBA_CMD="mamba"
        print_info "Using existing mamba installation"
    elif check_conda_installed; then
        # Check if conda has mamba installed
        if conda list mamba 2>/dev/null | grep -q "^mamba "; then
            MAMBA_CMD="mamba"
            print_info "Using mamba from existing conda installation"
        else
            print_warning "Conda found but mamba not installed"
            print_info "Installing mamba into base conda environment..."
            conda install -n base -c conda-forge mamba -y
            MAMBA_CMD="mamba"
        fi
    else
        print_warning "Mamba/Conda not found. Installing Miniforge..."
        install_mamba
        # Script will exit after installation, user needs to rerun
    fi
    
    # Check if environment file exists
    if [ ! -f "$ENV_YML" ]; then
        print_error "Environment file not found: $ENV_YML"
        exit 1
    fi
    
    print_info "Using environment file: $ENV_YML"
    
    # Check if environment already exists
    if $MAMBA_CMD env list | grep -q "^${ENV_NAME} "; then
        print_warning "Environment '$ENV_NAME' already exists"
        read -p "Do you want to remove and recreate it? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_step "Removing existing environment..."
            $MAMBA_CMD env remove -n "$ENV_NAME" -y
        else
            print_info "Updating existing environment..."
            $MAMBA_CMD env update -n "$ENV_NAME" -f "$ENV_YML"
            print_success "Environment updated successfully!"
            print_activation_instructions_mamba "$MAMBA_CMD"
            return 0
        fi
    fi
    
    # Create environment
    print_step "Creating mamba environment '$ENV_NAME'..."
    $MAMBA_CMD env create -f "$ENV_YML"
    
    print_success "Mamba environment created successfully!"
    print_activation_instructions_mamba "$MAMBA_CMD"
}

# Function to print mamba activation instructions
print_activation_instructions_mamba() {
    local mamba_cmd=$1
    echo ""
    print_info "To activate the environment, run:"
    echo -e "  ${CYAN}mamba activate $ENV_NAME${NC}"
    echo ""
    echo -e "  or"
    echo ""
    echo -e "  ${CYAN}conda activate $ENV_NAME${NC}"
}

# Function to setup venv environment
setup_venv_env() {
    print_step "Setting up Python virtual environment..."
    
    # Check if Python 3 is available
    if ! command -v python3 &> /dev/null; then
        print_error "Python 3 is not installed. Please install Python 3 first."
        exit 1
    fi
    
    PYTHON_VERSION=$(python3 --version | cut -d' ' -f2)
    print_info "Using Python $PYTHON_VERSION"
    
    # Check if requirements file exists
    if [ ! -f "$REQUIREMENTS_TXT" ]; then
        print_error "Requirements file not found: $REQUIREMENTS_TXT"
        exit 1
    fi
    
    print_info "Using requirements file: $REQUIREMENTS_TXT"
    
    # Environment directory
    VENV_DIR="$PROJECT_ROOT/$ENV_NAME"
    
    # Check if environment already exists
    if [ -d "$VENV_DIR" ]; then
        print_warning "Virtual environment already exists at: $VENV_DIR"
        read -p "Do you want to remove and recreate it? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_step "Removing existing virtual environment..."
            rm -rf "$VENV_DIR"
        else
            print_info "Using existing virtual environment"
            print_step "Upgrading pip..."
            "$VENV_DIR/bin/pip" install --upgrade pip
            
            print_step "Installing/updating packages from requirements..."
            "$VENV_DIR/bin/pip" install -r "$REQUIREMENTS_TXT"
            
            print_success "Virtual environment updated successfully!"
            print_activation_instructions_venv
            return 0
        fi
    fi
    
    # Create virtual environment
    print_step "Creating virtual environment at: $VENV_DIR"
    python3 -m venv "$VENV_DIR"
    
    # Upgrade pip
    print_step "Upgrading pip..."
    "$VENV_DIR/bin/pip" install --upgrade pip
    
    # Install packages
    print_step "Installing packages from requirements..."
    "$VENV_DIR/bin/pip" install -r "$REQUIREMENTS_TXT"
    
    print_success "Virtual environment created successfully!"
    print_activation_instructions_venv
}

# Function to print venv activation instructions
print_activation_instructions_venv() {
    echo ""
    print_info "To activate the environment, run:"
    echo -e "  ${CYAN}source $PROJECT_ROOT/$ENV_NAME/bin/activate${NC}"
    echo ""
    print_info "To deactivate, run:"
    echo -e "  ${CYAN}deactivate${NC}"
}

# Function to show help
show_help() {
    cat << EOF
${BLUE}COSMIC Python Environment Setup${NC}

${CYAN}Usage:${NC}
  $0 [--mamba | --venv]

${CYAN}Options:${NC}
  --mamba    Create environment using mamba/micromamba from cosmic-env.yml
  --venv     Create environment using Python venv from cosmic-requirements.txt
  --help     Show this help message

${CYAN}Examples:${NC}
  $0 --mamba
  $0 --venv

${CYAN}Description:${NC}
  This script sets up a Python environment for COSMIC development.
  
  ${YELLOW}Mamba mode:${NC}
    - Installs Miniforge (includes mamba) if not already installed
    - Creates a conda environment from cosmic-env.yml
    - Environment name: cosmic
  
  ${YELLOW}Venv mode:${NC}
    - Creates a Python virtual environment
    - Installs packages from cosmic-requirements.txt
    - Environment location: $PROJECT_ROOT/cosmic

EOF
}

# Main script
main() {
    # Parse arguments
    if [ $# -eq 0 ]; then
        print_error "No arguments provided"
        echo ""
        show_help
        exit 1
    fi
    
    case "$1" in
        --mamba)
            setup_mamba_env
            ;;
        --venv)
            setup_venv_env
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            echo ""
            show_help
            exit 1
            ;;
    esac
}

# Run main function
main "$@"
