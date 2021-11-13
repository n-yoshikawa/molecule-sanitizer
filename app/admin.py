from django.contrib import admin
from .models import Molecule, Project, Evaluation, Tag

# Register your models here.
admin.site.register(Molecule)
admin.site.register(Evaluation)
admin.site.register(Tag)
admin.site.register(Project)
