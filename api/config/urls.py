"""URL config for multi-dwpc API project."""

from django.urls import include, path

urlpatterns = [
    path("api/", include("dwpc_api.urls")),
    path("api-auth/", include("rest_framework.urls")),
]
