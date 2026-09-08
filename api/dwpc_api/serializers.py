from __future__ import annotations

from rest_framework import serializers

from src.multi_dwpc_query import DEFAULT_B, DEFAULT_PATH_Z_MIN


class BaseGeneTargetSerializer(serializers.Serializer):
    target_id = serializers.CharField(required=False, allow_blank=False)
    target_name = serializers.CharField(required=False, allow_blank=False)
    gene_ids = serializers.ListField(
        child=serializers.IntegerField(), required=False, allow_empty=False
    )
    gene_symbols = serializers.ListField(
        child=serializers.CharField(), required=False, allow_empty=False
    )
    genes_raw = serializers.CharField(required=False, allow_blank=False)

    def validate(self, attrs):
        has_target = bool(attrs.get("target_id") or attrs.get("target_name"))
        has_genes = bool(
            attrs.get("gene_ids") or attrs.get("gene_symbols") or attrs.get("genes_raw")
        )
        if not has_target:
            raise serializers.ValidationError("Provide target_id or target_name.")
        if not has_genes:
            raise serializers.ValidationError(
                "Provide one of gene_ids, gene_symbols, or genes_raw."
            )
        return attrs


class MetapathRankingRequestSerializer(BaseGeneTargetSerializer):
    b = serializers.IntegerField(required=False, min_value=2, default=DEFAULT_B)
    seed = serializers.IntegerField(required=False, default=42)


class MetapathDrilldownRequestSerializer(BaseGeneTargetSerializer):
    metapath = serializers.CharField()
    path_top_k = serializers.IntegerField(required=False, min_value=1, default=500)
    path_z_min = serializers.FloatField(required=False, default=DEFAULT_PATH_Z_MIN)
    debug = serializers.BooleanField(required=False, default=False)


class OverallSubgraphsRequestSerializer(BaseGeneTargetSerializer):
    b = serializers.IntegerField(required=False, min_value=2, default=DEFAULT_B)
    seed = serializers.IntegerField(required=False, default=42)
    pool_z_min = serializers.FloatField(required=False, default=1.65)
    top_n_paths = serializers.IntegerField(required=False, min_value=1, default=40)
    top_n_shared_per_hop = serializers.IntegerField(required=False, min_value=1, default=15)
    path_top_k = serializers.IntegerField(required=False, min_value=1, default=500)
    path_z_min = serializers.FloatField(required=False, default=DEFAULT_PATH_Z_MIN)


class ReferenceLookupSerializer(serializers.Serializer):
    q = serializers.CharField(required=False, allow_blank=True, default="")
    limit = serializers.IntegerField(required=False, min_value=1, max_value=5000, default=100)
