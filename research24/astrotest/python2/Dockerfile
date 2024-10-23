FROM continuumio/anaconda2

# Install necessary Python packages via conda
#RUN conda install -y numpy scipy matplotlib pandas

# Add this line to accept the suite changes in Debian repositories
#RUN apt-get update --allow-releaseinfo-change

# Update and install system packages (including GCC and GFortran)
#RUN apt-get install -y gcc gfortran && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy all files from the current directory to the container's working directory
COPY ./ ./

# Clean up to reduce image size
#RUN apt-get clean && rm -rf /var/lib/apt/lists/*

