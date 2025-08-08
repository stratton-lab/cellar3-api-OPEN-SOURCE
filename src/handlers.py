from django.http import JsonResponse, HttpResponseForbidden


def handler404json(request, exception):
    return JsonResponse({'error': 'This resource is not available'}, status=404)


def handler400json(request, exception: Exception):
    return JsonResponse({'error': 'Invalid user request', 'details': exception}, status=400)


def handler500json(request):
    return JsonResponse({'error': 'Server Error'}, status=500)


def forbidden_view(request):
    return HttpResponseForbidden("Access to the admin site is forbidden in the production environment.")
