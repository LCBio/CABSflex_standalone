import numpy as np
import numpy.typing as npt
from typing import Dict, List, Tuple, Union, Optional, Any
from typing_extensions import Literal
from CABS.trajectory import Trajectory


class Filter:
    def __init__(self, trajectory, num=1000):
        """
        Class for performing trajectory filtering according to chosen criteria (currently -- lowest energy).
        :param trajectory: trajectory.Trajectory instance to be clustered.
        :param num: int the number of models to be filtered out.
        """
        super().__init__()
        self.trajectory = trajectory
        self.N = num

    @staticmethod
    def mdl_fltr(mdls, enrgs, num=None):
        """
        Assisting method for filtering.
        :param mdls: list of models.
        :param enrgs: list of energies.
        :param num: int number of models to be filtered out.
        :return: list of indeces of the filtered models.
        """
        if num is None:
            num = len(mdls)
        low_energy_ndxs = np.argsort(enrgs)
        if len(mdls) <= num:
            filtered_ndx = low_energy_ndxs
        else:
            filtered_ndx = low_energy_ndxs[:num]
        return filtered_ndx

    def cabs_filter(self):
        """
        Default CABS-dock filtering method.
        :return: trajectory. Trajectory instance with filtered out models, list of indices of the models.
        """
        n_replicas = self.trajectory.coordinates.shape[0]
        n_models = self.trajectory.coordinates.shape[1]
        fromeach = int(self.N / n_replicas)
        filtered_models = []
        filtered_headers = []
        filtered_total_ndx = []

        for i, replica in enumerate(self.trajectory.coordinates):
            energies = [header.get_energy(number_of_peptides=self.trajectory.number_of_peptides) for header in
                        self.trajectory.headers if header.replica == i + 1]
            headers = [header for header in self.trajectory.headers if header.replica == i + 1]
            filtered_ndx = self.mdl_fltr(replica, energies, num=fromeach)
            if len(filtered_models) == 0:
                filtered_models = replica[filtered_ndx, :, :]
                filtered_headers = [headers[k] for k in filtered_ndx]
            else:
                filtered_models = np.concatenate([filtered_models, replica[filtered_ndx, :, :]])
                filtered_headers += [headers[k] for k in filtered_ndx]
            filtered_total_ndx.extend(np.array(filtered_ndx) + i * n_models)
        traj = Trajectory(self.trajectory.template, np.array([filtered_models]), filtered_headers)
        traj.number_of_peptides = self.trajectory.number_of_peptides
        return traj, filtered_total_ndx
