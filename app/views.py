from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.http import HttpResponse, JsonResponse
from django.shortcuts import redirect, render
from social_django.models import UserSocialAuth
from .models import UserProfile, Project, Molecule, Evaluation, Tag

import os
import threading
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
import pandas as pd
import urllib.parse

def create_new_user(user):
    profile = UserProfile.objects.create(user=user)
    Evaluation.objects.create(user=profile, evaluation_type='Like')
    Evaluation.objects.create(user=profile, evaluation_type='DisLike')
    return profile

def index(request):
    if request.user.is_authenticated:
        user = UserSocialAuth.objects.get(user_id=request.user.id)
        try:
            profile = UserProfile.objects.get(user=user)
        except UserProfile.DoesNotExist:
            profile = create_new_user(user)
        projects = Project.objects.filter(owner=profile)
        return render(request,'app/index.html', {'user': user, 'projects': projects})
    else:
        return render(request,'app/index.html')

@login_required
def new(request):
    if request.method == 'POST':
        smiles_list = request.POST['smiles_list']
        nlines = len(smiles_list.splitlines())
        if nlines > 10000:
            return render(request,'app/new.html', {'error': 'SMILES list is too long'})
        name = request.POST['project_name']
        user = UserSocialAuth.objects.get(user_id=request.user.id)
        profile = UserProfile.objects.get(user=user)
        project = Project.objects.create(owner=profile, raw_text=smiles_list, name=name)
        thread = threading.Thread(target=process_project, args=(project,))
        thread.start()
        return redirect('viewer', project_id=project.id)
    else:
        return render(request,'app/new.html')

def process_project(project):
    smiles_list = project.raw_text
    try:
        rule_of_five = Tag.objects.get(name='Ro5', badge_type='info')
    except:
        rule_of_five = Tag.objects.create(name='Ro5', badge_type='info')
    try:
        pains_tag = Tag.objects.get(name='PAINS', badge_type='danger')
    except:
        pains_tag = Tag.objects.create(name='PAINS', badge_type='danger')
    try:
        mcf_tag = Tag.objects.get(name='MCF', badge_type='danger')
    except:
        mcf_tag = Tag.objects.create(name='MCF', badge_type='danger')
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains_filter = FilterCatalog.FilterCatalog(params)
    _base_dir = os.path.split(__file__)[0]
    _mcf_csv = pd.read_csv(os.path.join(_base_dir, 'mcf.csv'))
    mcf_list = [Chem.MolFromSmarts(x) for x in _mcf_csv['smarts'].values]
    for smiles in smiles_list.splitlines():
        mol = Molecule.objects.create(project=project)
        mol.smiles = smiles
        mol_rdkit = Chem.MolFromSmiles(smiles)
        if mol is None:
            mol.save()
            continue
        inchikey = Chem.rdinchi.MolToInchiKey(mol_rdkit)
        mol.inchikey = inchikey
        path = os.path.join(os.path.dirname(__file__), f'static/app/mol_images/{mol.uuid}.png')
        d = rdMolDraw2D.MolDraw2DCairo(400, 300)
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol_rdkit)
        d.FinishDrawing()
        d.WriteDrawingText(path)
        # rule of five
        mol.molecular_weight = Descriptors.ExactMolWt(mol_rdkit)
        mol.logp = Descriptors.MolLogP(mol_rdkit)
        mol.h_bond_donor = Descriptors.NumHDonors(mol_rdkit)
        mol.h_bond_acceptors = Descriptors.NumHAcceptors(mol_rdkit)
        if mol.molecular_weight <= 500 and mol.logp <= 5 \
           and mol.h_bond_donor <= 5 and mol.h_bond_acceptors <= 10:
            mol.tags.add(rule_of_five)
        if pains_filter.HasMatch(mol_rdkit):
            mol.tags.add(pains_tag)
        h_mol = Chem.AddHs(mol_rdkit)
        if any(h_mol.HasSubstructMatch(smarts) for smarts in mcf_list):
            mol.tags.add(mcf_tag)
        mol.save()
    project.is_ready = True
    project.save()

def apply_filter(molecules, profile, filter_names, smiles):
    if not molecules:
        return molecules
    if filter_names:
        filter_names = filter_names.split('+')
        for filter_name in filter_names:
            if filter_name in ['Ro5', 'PAINS', 'MCF']:
                tag = Tag.objects.get(name=filter_name)
                molecules = molecules.filter(tags=tag)
            elif filter_name == 'like':
                evaluation = Evaluation.objects.get(user=profile, evaluation_type='Like')
                molecules = molecules.filter(evaluation=evaluation)
            elif filter_name == 'dislike':
                evaluation = Evaluation.objects.get(user=profile, evaluation_type='DisLike')
                molecules = molecules.filter(evaluation=evaluation)
    if smiles is not None:
        p = Chem.MolFromSmiles(smiles)
        if p:
            molecules = [m for m in molecules if Chem.MolFromSmiles(m.smiles).HasSubstructMatch(p)]
    return molecules


def viewer(request, project_id):
    project = Project.objects.get(id=project_id)
    n = 100
    if project.is_ready:
        molecules = Molecule.objects.filter(project=project)
        page = int(request.GET.get('page', default='1'))
        filter_names = request.GET.get('filter')
        smiles = request.GET.get('smiles')
        if request.user.is_authenticated:
            user = UserSocialAuth.objects.get(user_id=request.user.id)
            profile = UserProfile.objects.get(user=user)
        else:
            try:
                profile = UserProfile.objects.get(user=None)
            except:
                profile = UserProfile.objects.create(user=None)
        molecules = apply_filter(molecules, profile, filter_names, smiles)
        filter_names_url = None
        if filter_names:
            filter_names_url = filter_names.rstrip()
            filter_names = filter_names.split('+')
        pages = [i+1 for i in range((len(molecules) + n - 1) // n)]
        if len(pages) > 10:
            if page <= 5:
                pages = pages[:8] + [-1] + pages[-1:]
            elif page <= len(pages) - 5:
                pages = pages[:1] + [-1] + pages[page-4:page+3] + [-1] + pages[-1:]
            else:
                pages = pages[:1] + [-1] + pages[-8:]
        molecules = molecules[n*(page-1):n*page]
        mol_info = []
        for molecule in molecules:
            liked = True if len(molecule.evaluation.filter(user=profile, evaluation_type='Like')) > 0 else False
            disliked = True if len(molecule.evaluation.filter(user=profile, evaluation_type='DisLike')) > 0 else False
            mol_info.append((molecule, liked, disliked))
        url_params = None
        if filter_names_url and smiles:
            url_params = f'filter={urllib.parse.quote_plus(filter_names_url)}&smiles={urllib.parse.quote_plus(smiles)}'
        elif filter_names_url:
            url_params = f'filter={urllib.parse.quote_plus(filter_names_url)}'
        elif smiles:
            url_params = f'smiles={urllib.parse.quote_plus(smiles)}'
        return render(request,'app/viewer.html', 
                   {'project': project, 'mol_info': mol_info, 'project_id': project_id,
                    'filter_names': filter_names, 'filter_names_url': filter_names_url, 
                    'pages': pages, 'current_page': page, 'prev_page': page-1 if page != 1 else None,
                    'next_page': page+1 if pages and page != pages[-1] else None,
                    'smiles': smiles, 'url_params': url_params})
    else:
        return render(request,'app/viewer.html', {'project': project})


def export(request, project_id):
    project = Project.objects.get(id=project_id)
    if project.is_ready:
        molecules = Molecule.objects.filter(project=project)
        page = int(request.GET.get('page', default='1'))
        filter_names = request.GET.get('filter')
        smiles = request.GET.get('smiles')
        if request.user.is_authenticated:
            user = UserSocialAuth.objects.get(user_id=request.user.id)
            profile = UserProfile.objects.get(user=user)
        else:
            profile = UserProfile.objects.get(user=None)
        molecules = apply_filter(molecules, profile, filter_names, smiles)
        return HttpResponse("\n".join([mol.smiles for mol in molecules]))
    else:
        return HttpResponse("")

def evaluation(request):
    if request.method == 'POST' and request.is_ajax():
        evaluation_type = request.POST.get('evaluation_type')
        if request.user.is_authenticated:
            user = UserSocialAuth.objects.get(user_id=request.user.id)
            profile = UserProfile.objects.get(user=user)
        else:
            profile = UserProfile.objects.get(user=None)
        try:
            like = Evaluation.objects.get(user=profile, evaluation_type='Like')
        except:
            like = Evaluation.objects.create(user=profile, evaluation_type='Like')
        try:
            dislike = Evaluation.objects.get(user=profile, evaluation_type='DisLike')
        except:
            dislike = Evaluation.objects.create(user=profile, evaluation_type='DisLike')

        action = request.POST.get('action')
        if action == 'all':
            filter_names = request.POST.get('filter')
            if filter_names == "":
                filter_names = None
            smiles = request.POST.get('smiles')
            if smiles == "":
                smiles = None
            project_id = request.POST.get('project_id')
            project = Project.objects.get(id=project_id)
            molecules = Molecule.objects.filter(project=project)
            molecules = apply_filter(molecules, profile, filter_names, smiles)
            for molecule in molecules:
                with transaction.atomic():
                    try:
                        molecule.evaluation.remove(like)
                    except:
                        pass
                    try:
                        molecule.evaluation.remove(dislike)
                    except:
                        pass
                    if evaluation_type == 'like':
                        molecule.evaluation.add(like)
                    elif evaluation_type == 'dislike':
                        molecule.evaluation.add(dislike)
                    molecule.like_count = len(molecule.evaluation.filter(evaluation_type='Like'))
                    molecule.dislike_count = len(molecule.evaluation.filter(evaluation_type='DisLike'))
                    molecule.save()
            return JsonResponse({'success': False});

        molecule_id = request.POST.get('molecule_id')
        molecule = Molecule.objects.get(id=molecule_id)
        with transaction.atomic():
            if evaluation_type == 'like':
                if action == 'add':
                    molecule.evaluation.add(like)
                elif action == 'remove':
                    molecule.evaluation.remove(like)
                elif action == 'toggle':
                    molecule.evaluation.add(like)
                    molecule.evaluation.remove(dislike)
            elif evaluation_type == 'dislike':
                if action == 'add':
                    molecule.evaluation.add(dislike)
                elif action == 'remove':
                    molecule.evaluation.remove(dislike)
                elif action == 'toggle':
                    molecule.evaluation.add(dislike)
                    molecule.evaluation.remove(like)
            molecule.like_count = len(molecule.evaluation.filter(evaluation_type='Like'))
            molecule.dislike_count = len(molecule.evaluation.filter(evaluation_type='DisLike'))
            molecule.save()
        return JsonResponse({'like_count': molecule.like_count,
                             'dislike_count': molecule.dislike_count})
    else:
        return JsonResponse({'type': 'failure'})

def molecule(request, molecule_id):
    molecule = Molecule.objects.get(uuid=molecule_id)
    if molecule.tags.filter(name='PAINS'):
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        pains = FilterCatalog.FilterCatalog(params)
        mol_rdkit = Chem.MolFromSmiles(molecule.smiles)
        matches = pains.GetMatches(mol_rdkit)
        pains = []
        for i, match in enumerate(matches):
            hatoms = [x[1] for x in match.GetFilterMatches(mol_rdkit)[0].atomPairs]
            hbonds = []
            for bond in mol_rdkit.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if a1 in hatoms and a2 in hatoms:
                    hbonds.append(bond.GetIdx())
            path = os.path.join(os.path.dirname(__file__), f'static/app/mol_images/{molecule.uuid}-pains-{i}.png')
            d = rdMolDraw2D.MolDraw2DCairo(400, 300)
            rdMolDraw2D.PrepareAndDrawMolecule(d, mol_rdkit, highlightAtoms=hatoms, highlightBonds=hbonds)
            d.FinishDrawing()
            d.WriteDrawingText(path)
            pains.append(i)
        return render(request,'app/molecule.html', {'molecule': molecule, 'pains': pains})
    return render(request,'app/molecule.html', {'molecule': molecule})

def rename_project(request):
    if request.method == 'POST' and request.is_ajax():
        project_id = request.POST.get('project_id')
        new_name = request.POST.get('new_name')
        with transaction.atomic():
            project = Project.objects.get(id=project_id)
            project.name = new_name
            project.save()
        return JsonResponse({'success': True})
    else:
        return JsonResponse({'success': False})

def remove_project(request):
    if request.method == 'POST' and request.is_ajax():
        project_id = request.POST.get('project_id')
        with transaction.atomic():
            project = Project.objects.filter(id=project_id).delete()
        return JsonResponse({'success': True})
    else:
        return JsonResponse({'success': False})

def search(request, project_id):
    project = Project.objects.get(id=project_id)
    smiles = request.GET.get('smiles')
    if smiles is not None:
        smiles = smiles.rstrip()
    return render(request, 'app/editor.html', {'project': project, 'smiles': smiles})
