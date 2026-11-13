from __future__ import annotations

from collections.abc import Callable
from typing import Any

from drf_spectacular.utils import OpenApiParameter, extend_schema
from rest_framework import status, viewsets
from rest_framework.decorators import action
from rest_framework.exceptions import APIException, ValidationError
from rest_framework.response import Response
from rest_framework.serializers import Serializer

from src.null_bundle import BundleMismatchError

from .serializers import (
    BiologicalProcessListSerializer,
    BiologicalProcessSerializer,
    MetapathDrilldownRequestSerializer,
    MetapathDrilldownResponseSerializer,
    MetapathRankingResponseSerializer,
    OverallSubgraphsResponseSerializer,
    QueryRequestSerializer,
    RandomGenesSerializer,
    ReferenceLookupSerializer,
)
from .services import (
    list_biological_processes,
    random_biological_process,
    random_genes,
    run_metapath_drilldown,
    run_metapath_ranking,
    run_overall_subgraphs,
)


class NullBundleUnavailable(APIException):
    status_code = status.HTTP_503_SERVICE_UNAVAILABLE
    default_code = "null_bundle_unavailable"


def _run_query(
    serializer_class: type[Serializer],
    request,
    run: Callable[[dict[str, Any]], dict[str, Any]],
) -> Response:
    serializer = serializer_class(data=request.data)
    serializer.is_valid(raise_exception=True)
    try:
        payload = run(serializer.validated_data)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc
    except BundleMismatchError as exc:
        raise NullBundleUnavailable(str(exc)) from exc
    return Response(payload, status=status.HTTP_200_OK)


class QueryViewSet(viewsets.ViewSet):
    @extend_schema(request=QueryRequestSerializer, responses=MetapathRankingResponseSerializer)
    @action(detail=False, methods=["post"], url_path="metapath-ranking")
    def metapath_ranking(self, request):
        return _run_query(QueryRequestSerializer, request, run_metapath_ranking)

    @extend_schema(
        request=MetapathDrilldownRequestSerializer,
        responses=MetapathDrilldownResponseSerializer,
    )
    @action(detail=False, methods=["post"], url_path="metapath-drilldown")
    def metapath_drilldown(self, request):
        return _run_query(MetapathDrilldownRequestSerializer, request, run_metapath_drilldown)

    @extend_schema(request=QueryRequestSerializer, responses=OverallSubgraphsResponseSerializer)
    @action(detail=False, methods=["post"], url_path="overall-subgraphs")
    def overall_subgraphs(self, request):
        return _run_query(QueryRequestSerializer, request, run_overall_subgraphs)


class ReferenceViewSet(viewsets.ViewSet):
    @extend_schema(
        parameters=[
            OpenApiParameter("q", str, description="Substring of the name or GO identifier."),
            OpenApiParameter("limit", int),
        ],
        responses=BiologicalProcessListSerializer,
    )
    @action(detail=False, methods=["get"], url_path="biological-processes")
    def biological_processes(self, request):
        serializer = ReferenceLookupSerializer(data=request.query_params)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        results = list_biological_processes(q=data["q"], limit=data["limit"])
        return Response({"results": results}, status=status.HTTP_200_OK)

    @extend_schema(responses=BiologicalProcessSerializer)
    @action(detail=False, methods=["get"], url_path="biological-processes/random")
    def random_biological_process(self, request):
        return Response(random_biological_process(), status=status.HTTP_200_OK)

    @extend_schema(responses=RandomGenesSerializer)
    @action(detail=False, methods=["get"], url_path="genes/random")
    def random_genes(self, request):
        return Response({"genes": random_genes()}, status=status.HTTP_200_OK)
