
# save path
import pickle
import os

def SaveObj(obj, filename = None):
    """Save object into pkl-formatted pickle file.

    Parameters
    ----------
    obj : :obj:`~coolname.Budoids_class.Budoid`
        Budoid object
    filename : :class:`str`
        File name of data file.

    """

    if filename is None:
        filename = os.path.join(obj.outdir, 'Budoids.pkl')
    with open(filename, "wb") as fn:
        pickle.dump(obj, fn)


def ReadObj(filename):
    r"""Read .pkl-formatted pickle file.

    Parameters
    ----------
    filename : :class:`str`
        File name of data file.

    """

    with open(filename, "rb") as fn:
        return pickle.load(fn)


# def SavePath(path, out_dir = '.'):
#     fn = os.path.join(out_dir, 'path.pkl')

#     with open(fn, "wb") as fp:
#         pickle.dump(path, fp)
