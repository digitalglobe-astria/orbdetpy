#!/bin/bash
# gen_rpc_stubs.sh - Generate Python gRPC stubs.
# Copyright (C) 2019 University of Texas
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

BASEDIR=$(dirname "$0")
pushd $BASEDIR > /dev/null

python -m grpc_tools.protoc -I src/main/proto/ --python_out=./protobuf --grpc_python_out=./protobuf src/main/proto/*.proto

# Ugly fix for module import path issue in gRPC generated Python stub
sed -i "s/import /import orbdetpy.protobuf./g" ./protobuf/*_grpc.py
sed -i "s/import orbdetpy.protobuf.grpc/import grpc/g" ./protobuf/*_grpc.py

popd > /dev/null
