from __future__ import annotations

from rest_framework import serializers


class QueryRequestSerializer(serializers.Serializer):
    target_id = serializers.CharField()
    genes_raw = serializers.CharField(
        help_text="Gene symbols or Entrez IDs, one per line or comma-separated."
    )


class MetapathDrilldownRequestSerializer(QueryRequestSerializer):
    metapath = serializers.CharField()


class ReferenceLookupSerializer(serializers.Serializer):
    q = serializers.CharField(required=False, allow_blank=True, default="")
    limit = serializers.IntegerField(required=False, min_value=1, max_value=5000, default=100)


class BiologicalProcessSerializer(serializers.Serializer):
    identifier = serializers.CharField()
    name = serializers.CharField()


class BiologicalProcessListSerializer(serializers.Serializer):
    results = BiologicalProcessSerializer(many=True)


class RandomGenesSerializer(serializers.Serializer):
    genes = serializers.ListField(child=serializers.CharField())


class MetapathRankingRowSerializer(serializers.Serializer):
    metapath = serializers.CharField()
    real_mean_score = serializers.FloatField()
    null_mean_score = serializers.FloatField()
    null_std_score = serializers.FloatField(allow_null=True)
    diff = serializers.FloatField()
    effect_size_z = serializers.FloatField(allow_null=True)
    p_value = serializers.FloatField(allow_null=True)


class MetapathRankingResponseSerializer(serializers.Serializer):
    target_id = serializers.CharField()
    gene_ids = serializers.ListField(child=serializers.IntegerField())
    metapath_ranking = MetapathRankingRowSerializer(many=True)
    elapsed_s = serializers.FloatField()


class IntermediateSharingRowSerializer(serializers.Serializer):
    intermediate_id = serializers.CharField()
    intermediate = serializers.CharField()
    n_genes_sharing = serializers.IntegerField()
    gene_ids = serializers.CharField(help_text="Comma-separated Entrez IDs.")
    genes = serializers.CharField(help_text="Comma-separated gene symbols.")


class HeatmapAxisItemSerializer(serializers.Serializer):
    id = serializers.CharField()
    name = serializers.CharField()


class HeatmapSerializer(serializers.Serializer):
    intermediates = HeatmapAxisItemSerializer(many=True)
    genes = HeatmapAxisItemSerializer(many=True)
    cells = serializers.ListField(
        child=serializers.ListField(child=serializers.IntegerField()),
        help_text="Rows follow intermediates, columns follow genes; 1 if the gene uses it.",
    )


class SubgraphNodeSerializer(serializers.Serializer):
    id = serializers.CharField()
    node_id = serializers.CharField()
    name = serializers.CharField()
    type = serializers.CharField()
    column = serializers.IntegerField()
    label = serializers.CharField()
    scale = serializers.FloatField()


class SubgraphEdgeSerializer(serializers.Serializer):
    source = serializers.CharField()
    target = serializers.CharField()


class SubgraphSerializer(serializers.Serializer):
    title = serializers.CharField()
    columns = serializers.ListField(child=serializers.CharField())
    nodes = SubgraphNodeSerializer(many=True)
    edges = SubgraphEdgeSerializer(many=True)


class MetapathDrilldownResponseSerializer(serializers.Serializer):
    target_id = serializers.CharField()
    gene_ids = serializers.ListField(child=serializers.IntegerField())
    metapath = serializers.CharField()
    path_z_min = serializers.FloatField()
    intermediate_sharing = IntermediateSharingRowSerializer(many=True)
    heatmap = HeatmapSerializer(allow_null=True)
    paths = serializers.ListField(
        child=serializers.DictField(),
        help_text="gene, gene_id, metapath, path_rank, path_score and hop_N_id/type/name.",
    )
    subgraph = SubgraphSerializer(allow_null=True)
    diagnostics = serializers.DictField()


class OverallSubgraphsResponseSerializer(serializers.Serializer):
    target_id = serializers.CharField()
    gene_ids = serializers.ListField(child=serializers.IntegerField())
    pooled_metapaths = serializers.ListField(child=serializers.CharField())
    pool_z_min = serializers.FloatField()
    pool_z_used = serializers.FloatField()
    top_shared_subgraph = SubgraphSerializer(allow_null=True)
    top_paths_subgraph = SubgraphSerializer(allow_null=True)
