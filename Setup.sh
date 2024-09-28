#!/bin/bash

proj_dir="$(pwd)"
my_code="$(pwd)/code"
my_data="$(pwd)/data"
my_results="$(pwd)/results"
my_softwares="$(pwd)/softwares"
my_docs="$(pwd)/docs"
mkdir -p $my_code
mkdir -p $my_data
mkdir -p $my_results
mkdir -p $my_softwares
mkdir -p $my_docs
env_file="$my_code/project_env.sh"

cat <<EOL > "$env_file"
#!/bin/bash
export proj_dir="$(pwd)"
export my_code="$(pwd)/code"
export my_data="$(pwd)/data"
export my_results="$(pwd)/results"
export my_softwares="$(pwd)/softwares"
export my_docs="$(pwd)/docs"
export LAMPLD_dir="$(pwd)/softwares/LAMPLD_perl/"
EOL

# Make the file executable
chmod +x "$env_file"


cd $my_code

