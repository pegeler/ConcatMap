from importlib import metadata

from concatmap.mapper import concatmap

__version__ = metadata.version('concatmap')
__all__ = ['concatmap']
