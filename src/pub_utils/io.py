# import h5py
# import zarr
import pickle


def handle_pickle(data=None, filename="data.pkl", mode="save"):
    if mode == "save":
        if data is None:
            raise ValueError("You must provide data to save.")
        with open(filename, "wb") as f:
            pickle.dump(data, f)
        print(f"Successfully saved to {filename}")
        
    elif mode == "load":
        try:
            with open(filename, "rb") as f:
                return pickle.load(f)
        except FileNotFoundError:
            print(f"Error: The file {filename} was not found.")
            return None
    else:
        print("Invalid mode. Please use 'save' or 'load'.")