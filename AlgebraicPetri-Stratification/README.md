# GTRI Model Stratification using AlgebraicPetri

## Set-up Environment

1. Have Julia downloaded and installed
1. Clone repository
1. From this directory run: `julia --project setup.jl`
1. Run the example script: `julia --project stratification-notebook.jl`

## Included Files

- `setup.jl`: A simple Julia script to help set up a local environment
- `models/`: A collection of example `json` formatted Petri Net models
- `connection_graphs/`: A collection of example `json` formatted interconnection
  graphs for defining stratified models
- `ModelStratify.jl`: Julia module that defines the API methods to interact with
  `AlgebraicPetri` and `Catlab` libraries
- `stratification-notebook.jl`: An example script/notebook of using the API to
  stratify models either demographically or spatially

## `ModelStratify.jl` API Methods

- `diff_strat(epi_model, connection_graph)`

    This function takes in a `LabelledPetriNet` and a graph which describes
    geographical connections. It returns a `LabelledPetriNet` which models
    diffusion between geographic populations described by the given graph.

- `dem_strat(epi_model, connection_graph, sus_state, exp_state, inf_states)`

    This function takes in a `LabelledPetriNet` and a graph which describes
    infection connections between populations. It also takes in the symbol used
    for susceptible states, the symbol used for the exposed state, and an array
    of symbols for states that can expose susceptible states. It returns a
    `LabelledPetriNet` which models diffusion between populations described by the
    given graph.

- `serialize(x; io)`

    Serialize an ACSet object to a JSON string

- `deserialize(x, ACSetType)`

    Deserialize a JSON string to an object of the given ACSet type

- `save_petri(epi_model, file_name, image_format)`

    This function takes an `AlgebraicPetri` model and saves it as an image given
    a filename and an image format (e.g. "png", "svg", etc.)

- `save_json(epi_model, file_name)`

    This function takes an `AlgebraicPetri` model and saves it as a `.json` file
    given a filename

- `save_model(epi_model, file_name)`

    This function is a shortcut to generate both an SVG file of the Petri model
    graph as well as the json file of the serialization of the ACSet object
    given an `AlgebraicPetri` model and a filename.
