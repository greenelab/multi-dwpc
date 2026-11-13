from rest_framework.routers import SimpleRouter

from .views import QueryViewSet, ReferenceViewSet

router = SimpleRouter()
router.register("queries", QueryViewSet, basename="queries")
router.register("reference", ReferenceViewSet, basename="reference")

urlpatterns = router.urls
