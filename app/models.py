from django.db import models
from social_django.models import UserSocialAuth

import uuid


class UserProfile(models.Model):
    user = models.OneToOneField(UserSocialAuth, on_delete=models.CASCADE, null=True)


class Project(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4)
    owner = models.ForeignKey(UserProfile, on_delete=models.CASCADE)
    name = models.CharField(max_length=64)
    is_ready = models.BooleanField(default=False)
    raw_text = models.TextField(default="")
    created_at = models.DateField(auto_now_add=True)
    updated_at = models.DateField(auto_now_add=True)


class Evaluation(models.Model):
    user = models.ForeignKey(UserProfile, on_delete=models.CASCADE)
    LIKE = 'LIKE'
    DISLIKE = 'DISLIKE'
    REACTION_CHOICES = [(LIKE, 'Like'), (DISLIKE, 'Dislike')]
    evaluation_type = models.CharField(
        max_length = 16,
        choices = REACTION_CHOICES,
        default = LIKE,
    )

class EvaluationDetail(models.Model):
    user = models.ForeignKey(UserProfile, on_delete=models.CASCADE)
    score1 = models.IntegerField(default=-1)
    score2 = models.IntegerField(default=-1)
    score3 = models.IntegerField(default=-1)
    score4 = models.IntegerField(default=-1)
    score5 = models.IntegerField(default=-1)
    free_form = models.TextField(default="")

class Tag(models.Model):
    name = models.CharField(max_length=32)

class Molecule(models.Model):
    uuid = models.UUIDField(default=uuid.uuid4)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    smiles = models.TextField()
    inchikey = models.CharField(max_length=27)
    evaluation = models.ManyToManyField(Evaluation)
    evaluation_detail = models.ManyToManyField(EvaluationDetail)
    tags = models.ManyToManyField(Tag)
    like_count = models.IntegerField(default=0)
    dislike_count = models.IntegerField(default=0)
    molecular_weight = models.FloatField(default=0)
    logp = models.FloatField(default=0)
    h_bond_donor = models.IntegerField(default=0)
    h_bond_acceptors = models.IntegerField(default=0)
