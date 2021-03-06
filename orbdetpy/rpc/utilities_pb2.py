# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: utilities.proto

import sys
_b=sys.version_info[0]<3 and (lambda x:x) or (lambda x:x.encode('latin1'))
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


import orbdetpy.rpc.messages_pb2 as messages__pb2


DESCRIPTOR = _descriptor.FileDescriptor(
  name='utilities.proto',
  package='',
  syntax='proto3',
  serialized_options=_b('\n\016org.astria.rpcB\020UtilitiesRequestP\000'),
  serialized_pb=_b('\n\x0futilities.proto\x1a\x0emessages.proto2@\n\tUtilities\x12\x33\n\timportTDM\x12\x0f.ImportTDMInput\x1a\x13.Measurement2DArray\"\x00\x42$\n\x0eorg.astria.rpcB\x10UtilitiesRequestP\x00\x62\x06proto3')
  ,
  dependencies=[messages__pb2.DESCRIPTOR,])



_sym_db.RegisterFileDescriptor(DESCRIPTOR)


DESCRIPTOR._options = None

_UTILITIES = _descriptor.ServiceDescriptor(
  name='Utilities',
  full_name='Utilities',
  file=DESCRIPTOR,
  index=0,
  serialized_options=None,
  serialized_start=35,
  serialized_end=99,
  methods=[
  _descriptor.MethodDescriptor(
    name='importTDM',
    full_name='Utilities.importTDM',
    index=0,
    containing_service=None,
    input_type=messages__pb2._IMPORTTDMINPUT,
    output_type=messages__pb2._MEASUREMENT2DARRAY,
    serialized_options=None,
  ),
])
_sym_db.RegisterServiceDescriptor(_UTILITIES)

DESCRIPTOR.services_by_name['Utilities'] = _UTILITIES

# @@protoc_insertion_point(module_scope)
