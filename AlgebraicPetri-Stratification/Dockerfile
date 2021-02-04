FROM julia

RUN apt-get update && apt-get install -y \
	gcc \
	graphviz \
	g++ \
	make \
	python3-pip

RUN pip3 install jupyter

COPY setup.jl .
COPY Project.toml .
RUN julia --project setup.jl

COPY build_sysimage.jl .
RUN julia --project build_sysimage.jl

COPY . .

# EXPOSE 8889
# ENTRYPOINT jupyter notebook --no-browser --allow-root --ip=0.0.0.0 --port=8889

ENTRYPOINT [ "julia", "-JSysImage.so", "--project", "cli.jl" ]
