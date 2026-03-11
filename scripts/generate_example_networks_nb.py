import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    """
    Generate and visualize example networks that showcase each network metric.

    Saves 2-panel PNGs (adjacency heatmap + network diagram) for each of the
    5 metrics defined in pub_utils.metrics:
        1. density   — sparse Erdos-Renyi graph
        2. modularity — stochastic block model with clear communities
        3. hubness   — star-like hub-dominated graph
        4. reciprocity — symmetrized graph with reciprocal edges
        5. clustering — triangle-rich graph with high transitivity

    Reads:  nothing (generates synthetic graphs)
    Saves:  connectomes/metrics/example_sparse.png
            connectomes/metrics/example_modular.png
            connectomes/metrics/example_hub.png
            connectomes/metrics/example_reciprocal.png
            connectomes/metrics/example_clustered.png

    Usage:
        python scripts/generate_example_networks.py
    """

    import matplotlib
    matplotlib.use("Agg")

    from pub_utils.generate import (
        generate_sparse,
        generate_modular,
        generate_hub,
        generate_reciprocal,
        generate_clustered,
    )

    OUT_DIR = "connectomes/metrics"

    GENERATORS = [
        ("sparse",     generate_sparse,     {"N": 8, "E": 10}),
        ("modular",    generate_modular,    {"N": 10, "E": 20}),
        ("hub",        generate_hub,        {"N": 8, "E": 14}),
        ("reciprocal", generate_reciprocal, {"N": 8, "E": 14}),
        ("clustered",  generate_clustered,  {"N": 8, "E": 18}),
    ]

    for name, fn, kwargs in GENERATORS:
        save_path = f"{OUT_DIR}/example_{name}.png"
        result = fn(save_path=save_path, **kwargs)
        print(f"  {name:12s}  {result['metric_name']} = {result['metric_value']:.3f}  -> {save_path}")

    print("\nDone.")

    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
