# Generated by the gRPC Python protocol compiler plugin. DO NOT EDIT!
import grpc

import orbdetpy.rpc.messages_pb2 as messages__pb2


class SimulationStub(object):
  # missing associated documentation comment in .proto file
  pass

  def __init__(self, channel):
    """Constructor.

    Args:
      channel: A grpc.Channel.
    """
    self.simulateMeasurements = channel.unary_unary(
        '/Simulation/simulateMeasurements',
        request_serializer=messages__pb2.Settings.SerializeToString,
        response_deserializer=messages__pb2.MeasurementArray.FromString,
        )


class SimulationServicer(object):
  # missing associated documentation comment in .proto file
  pass

  def simulateMeasurements(self, request, context):
    # missing associated documentation comment in .proto file
    pass
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')


def add_SimulationServicer_to_server(servicer, server):
  rpc_method_handlers = {
      'simulateMeasurements': grpc.unary_unary_rpc_method_handler(
          servicer.simulateMeasurements,
          request_deserializer=messages__pb2.Settings.FromString,
          response_serializer=messages__pb2.MeasurementArray.SerializeToString,
      ),
  }
  generic_handler = grpc.method_handlers_generic_handler(
      'Simulation', rpc_method_handlers)
  server.add_generic_rpc_handlers((generic_handler,))
