# Generated by the gRPC Python protocol compiler plugin. DO NOT EDIT!
import grpc

import orbdetpy.rpc.messages_pb2 as messages__pb2


class EstimationStub(object):
  # missing associated documentation comment in .proto file
  pass

  def __init__(self, channel):
    """Constructor.

    Args:
      channel: A grpc.Channel.
    """
    self.determineOrbit = channel.unary_unary(
        '/Estimation/determineOrbit',
        request_serializer=messages__pb2.DetermineOrbitInput.SerializeToString,
        response_deserializer=messages__pb2.EstimationOutputArray.FromString,
        )


class EstimationServicer(object):
  # missing associated documentation comment in .proto file
  pass

  def determineOrbit(self, request, context):
    # missing associated documentation comment in .proto file
    pass
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')


def add_EstimationServicer_to_server(servicer, server):
  rpc_method_handlers = {
      'determineOrbit': grpc.unary_unary_rpc_method_handler(
          servicer.determineOrbit,
          request_deserializer=messages__pb2.DetermineOrbitInput.FromString,
          response_serializer=messages__pb2.EstimationOutputArray.SerializeToString,
      ),
  }
  generic_handler = grpc.method_handlers_generic_handler(
      'Estimation', rpc_method_handlers)
  server.add_generic_rpc_handlers((generic_handler,))
