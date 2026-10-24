"""Semantic signature for retaining approved Stage 3 empty-transport caches.

This checks the independent sampling/count kernels, protocol primitives, study
inputs, ordered channels and timed circuits, rather than accepting an arbitrary
changed source hash. Transport-only extensions must not invalidate old data.
"""

import ast
from dataclasses import asdict, replace
import hashlib
import inspect
import json


def _code(fn, name=None):
    import textwrap
    node = ast.parse(textwrap.dedent(inspect.getsource(fn)))
    if name is not None:
        node.body[0].name = name
    return ast.dump(node, include_attributes=False)


def sampling_manifest():
    import ion_params as params
    import ion_timing as timing
    import steane_ft as ft
    import steane_noise as noise
    from ion_noise import IdleIonNoise, IonNoise
    code = {}
    for module, names in [(ft, ('stabilizer_circuit', 'run_ops', 'fault_list', 'ec_round',
                               'flagged_decode', 'standard_decode', 'syndrome_value')),
                          (noise, ('run_ops_frame', 'ec_round_frame', '_apply_payload', 'logical_failure',
                                   'physical_variants', 'all_locations', 'loc_type', 'validated_channel',
                                   'trigger_contributors', '_parity')),
                          (timing, ('idle_probability', 'operation_duration', 'insert_idles', 'timed_circuits'))]:
        for name in names:
            code[module.__name__ + '.' + name] = _code(getattr(module, name))
    for name in ('monte_carlo', 'count_faults'):
        kernel = getattr(noise, '_' + name + '_independent', getattr(noise, name))
        code[name] = _code(kernel, name)
        if kernel is not getattr(noise, name):
            # Exact opt-in dispatcher expected from Stage 4. A future change to
            # the default path must invalidate the compatibility signature.
            if name == 'count_faults':
                expected = '''def count_faults(noise=None):
    if getattr(noise, 'dynamic_channels', False):
        from ion_transport import conditional_count_faults
        return conditional_count_faults(noise)
    return _count_faults_independent(noise)
'''
            else:
                expected = '''def monte_carlo(noise, shots, rng, *, return_details=False):
    if getattr(noise, 'dynamic_channels', False):
        from ion_transport import conditional_monte_carlo
        return conditional_monte_carlo(noise, shots, rng, return_details=return_details)
    return _monte_carlo_independent(noise, shots, rng, return_details=return_details)
'''
            actual = _code(getattr(noise, name))
            if actual != ast.dump(ast.parse(expected), include_attributes=False):
                code[name + '_changed_dispatch'] = actual
    models = [IonNoise(p) for p in params.PRESETS.values()]
    j = params.JONES_MURALI_2026
    models += [IdleIonNoise(replace(j, T2=t)) for t in params.T2_SWEEP]
    models += [IdleIonNoise(replace(j, measurement_duration=t)) for t in params.MEASUREMENT_DURATION_SWEEP]
    models += [noise.UniformNoise(p) for p in params.UNIFORM_CALIBRATION_P]
    tables = []
    for model in models:
        circuits = getattr(model, 'circuits', noise.CIRCUITS)
        locs = noise.all_locations(circuits)
        tables.append({'locations': locs, 'channels': [list(model.channel(loc).items()) for loc in locs],
                       'dynamic': getattr(model, 'dynamic_channels', False)})
    data = {'code': code, 'models': tables, 'presets': [asdict(p) for p in params.PRESETS.values()],
            'round': ft.ROUND, 'supports': ft.SUPPORTS, 'data': ft.DATA,
            'trigger_types': noise.TRIGGER_TYPES,
            'study': [params.T2_SWEEP, params.MEASUREMENT_DURATION_SWEEP, params.UNIFORM_CALIBRATION_P]}
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
