FROM julia:1.5

RUN apt-get update && apt-get install -y \
 	gcc

COPY setup.jl .
COPY Project.toml .
RUN julia --project setup.jl

COPY build_sysimage.jl .
RUN julia --project build_sysimage.jl

COPY . .

EXPOSE 8001
ENTRYPOINT [ "julia", "-JSysImage.so", "--project", "server.jl" ]
