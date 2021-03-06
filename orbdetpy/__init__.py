# __init__.py - orbdetpy package initialization.
# Copyright (C) 2018-2019 University of Texas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import io
import os
import glob
import json
from .version import __version__
from orbdetpy.rpc.server import RemoteServer

def read_param(param):
    if (isinstance(param, str)):
        if (os.path.isfile(param)):
            with open(param, "r") as fp:
                data = json.load(fp)
        else:
            data = json.loads(param)
    elif (isinstance(param, io.TextIOBase)):
        data = json.load(param)
    else:
        data = param

    return(data)

def write_output_file(outfile, data):
    if (isinstance(outfile, str)):
        with open(outfile, "w") as fp:
            json.dump(data, fp)
    elif (isinstance(outfile, io.TextIOBase)):
        json.dump(data, outfile)

if (__name__ != '__main__'):
    _rootdir = os.path.dirname(os.path.abspath(__file__))
    _datadir = os.path.join(_rootdir, "data")
    _libsdir = os.path.join(_rootdir, "target")
    _jarfile = glob.glob(os.path.join(_libsdir, "astria*.jar"))[0]
    RemoteServer.connect(_datadir, _jarfile)
