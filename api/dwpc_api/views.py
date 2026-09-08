from __future__ import annotations

from collections import OrderedDict

from rest_framework import status, viewsets
from rest_framework.views import APIView
from rest_framework.decorators import action
from rest_framework.exceptions import ValidationError
from rest_framework.reverse import reverse
from rest_framework.response import Response

from .serializers import (
    MetapathDrilldownRequestSerializer,
    MetapathRankingRequestSerializer,
    OverallSubgraphsRequestSerializer,
    ReferenceLookupSerializer,
)
from .services import (
    list_biological_processes,
    list_genes,
    run_metapath_drilldown,
    run_metapath_ranking,
    run_overall_subgraphs,
)


class ApiRootView(APIView):
    """Browsable API landing page with concrete endpoint links."""

    def get(self, request):
        return Response(
            OrderedDict(
                [
                    ("metapath_ranking", reverse("queries-metapath-ranking", request=request)),
                    ("metapath_drilldown", reverse("queries-metapath-drilldown", request=request)),
                    ("overall_subgraphs", reverse("queries-overall-subgraphs", request=request)),
                    (
                        "reference_biological_processes",
                        reverse("reference-biological-processes", request=request),
                    ),
                    ("reference_genes", reverse("reference-genes", request=request)),
                ]
            )
        )


class QueryViewSet(viewsets.ViewSet):
    @action(detail=False, methods=["post"], url_path="metapath-ranking")
    def metapath_ranking(self, request):
        serializer = MetapathRankingRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        try:
            payload = run_metapath_ranking(serializer.validated_data)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        return Response(payload, status=status.HTTP_200_OK)

    @action(detail=False, methods=["post"], url_path="metapath-drilldown")
    def metapath_drilldown(self, request):
        serializer = MetapathDrilldownRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        try:
            payload = run_metapath_drilldown(serializer.validated_data)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        return Response(payload, status=status.HTTP_200_OK)

    @action(detail=False, methods=["post"], url_path="overall-subgraphs")
    def overall_subgraphs(self, request):
        serializer = OverallSubgraphsRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        try:
            payload = run_overall_subgraphs(serializer.validated_data)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        return Response(payload, status=status.HTTP_200_OK)


class ReferenceViewSet(viewsets.ViewSet):
    @action(detail=False, methods=["get"], url_path="biological-processes")
    def biological_processes(self, request):
        serializer = ReferenceLookupSerializer(data=request.query_params)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        payload = {
            "results": list_biological_processes(
                q=data.get("q", ""),
                limit=data.get("limit", 100),
            )
        }
        return Response(payload, status=status.HTTP_200_OK)

    @action(detail=False, methods=["get"], url_path="genes")
    def genes(self, request):
        serializer = ReferenceLookupSerializer(data=request.query_params)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        payload = {
            "results": list_genes(
                q=data.get("q", ""),
                limit=data.get("limit", 100),
            )
        }
        return Response(payload, status=status.HTTP_200_OK)
