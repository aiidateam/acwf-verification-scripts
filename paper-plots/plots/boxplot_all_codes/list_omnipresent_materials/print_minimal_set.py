# %%
import json
import pathlib as pl
import os
import sys

def get_list(set_names):
    """
    """
    DATA_FOLDER = "../../../code-data"
    with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
        labels_data = json.load(fhandle)
    
    reference_data_files = labels_data['references']['all-electron average']
    code_labels = labels_data['methods-main'].keys()

    for SET_NAME in set_names:
        try:
            with open(os.path.join(DATA_FOLDER, reference_data_files[SET_NAME])) as fhandle:
                ref_plugin_data = json.load(fhandle)
        except OSError as exc:
            raise ValueError(f"Error loading reference data found (set '{SET_NAME}')")
        ref_BM_fit_data = ref_plugin_data['BM_fit_data']
        ref_systems = set(
            el_conf for el_conf in ref_BM_fit_data.keys() if ref_BM_fit_data[el_conf] is not None
        )
        common_data = ref_systems
        print(f'{SET_NAME}', len(ref_systems))

        for code_label in code_labels:
            try:
                with open(os.path.join(DATA_FOLDER, labels_data['methods-main'][code_label][SET_NAME])) as fhandle:
                    plugin_data = json.load(fhandle)
            except OSError:
                print(f"No data found for {code_label} (set '{SET_NAME}')")
                sys.exit(1)

            plugin_BM_fit_data = plugin_data['BM_fit_data']

            plugin_systems = set(
                el_conf for el_conf in plugin_BM_fit_data.keys()
                if plugin_BM_fit_data[el_conf] is not None
            )
            # Take the systems that are both in the reference and plugin sets
            common_data = common_data.intersection(plugin_systems)
            print(f'{SET_NAME}  {code_label}', len(plugin_systems))

        print(common_data)

if __name__ == "__main__":
    get_list(['unaries', 'oxides'])

