from django.urls import path
from rest_framework.routers import DefaultRouter

from .views import ApiRootView, QueryViewSet, ReferenceViewSet

router = DefaultRouter()
router.include_root_view = False
# router.include_root_view = True
router.register("queries", QueryViewSet, basename="queries")
router.register("reference", ReferenceViewSet, basename="reference")

urlpatterns = [
	path("", ApiRootView.as_view(), name="api-root"),
]
urlpatterns += router.urls
