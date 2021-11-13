from django.urls import path, include
import django.contrib.auth.views

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('new/', views.new, name='new'),
    path('search/<uuid:project_id>', views.search, name='search'),
    path('viewer/<uuid:project_id>', views.viewer, name='viewer'),
    path('export/<uuid:project_id>', views.export, name='export'),
    path('molecule/<uuid:molecule_id>', views.molecule, name='molecule'),
    path('evaluation/', views.evaluation, name='evaluation'),
    path('rename_project', views.rename_project, name='rename_project'),
    path('remove_project', views.remove_project, name='remove_project'),
    path('logout/',
        django.contrib.auth.views.LogoutView.as_view(template_name = 'app/logout.html'),
        name='logout'),
]
