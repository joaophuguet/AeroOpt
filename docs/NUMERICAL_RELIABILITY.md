# Diagnósticos de confiabilidade numérica — v0.3.2

`aeroopt.verification` recebe histórias explícitas de refinamento. Não executa
solvers, não altera `SolverResult` e não fornece um limite matemático de erro.
Cada coeficiente é analisado separadamente; ordenar CM significa comparar seus
valores algébricos, sem atribuir preferência de projeto ao maior valor.

## Dados e uso

```python
from aeroopt.verification import RefinementRecord, RefinementStudy, compare_pair

def candidate(name, offset):
    records = tuple(
        RefinementRecord(method, n, cl + offset, None)
        for method, values in (
            ("hess_smith", (1.10, 1.04, 1.02)),
            ("linear_vortex", (1.08, 1.03, 1.015)),
        )
        for n, cl in zip((40, 80, 160), values, strict=True)
    )
    return RefinementStudy(name, records)

a, b = candidate("a", 0), candidate("b", .20)
history = a.analyze("cl")
comparison = compare_pair(a, b, "cl")
assert comparison["classification"] == "numerically_distinguishable"
assert compare_pair(a, b, "cm")["classification"] == "numerically_unresolved"
```

As quantidades de painéis são totais. São exigidos pelo menos três níveis,
estritamente crescentes por método, iguais entre métodos e candidatos. Não
há ordenação, interpolação nem preenchimento automático de histórias. Valores
ausentes ou registros `valid=False` permanecem `None`. A validade deve refletir
os gates do solver que gerou o registro. Rejeitam-se coeficientes não finitos,
níveis duplicados/decrescentes e políticas inválidas.

Limitações devem ser transmitidas explicitamente pelo chamador, por exemplo
`RefinementStudy(name, records, limitations=("FINITE_TE_UNVERIFIED",))`.
Essa API não inspeciona geometrias nem descobre automaticamente flags de um
resultado externo. Limitações explícitas impedem uma classificação resolvida.

## Critérios implementados

Para um coeficiente q e candidato c, denote por q(m,k) o valor do método m
no nível k. A história conserva cada incremento
`d(m,k) = q(m,k) - q(m,k-1)`.
O indicador é:

```text
U(c) = max(max_m |d(m,K)|, max_m q(m,K) - min_m q(m,K))
```

Comparam-se os dois últimos incrementos de cada método. Crescimento de módulo
acima de `|d_anterior|*(1+1e-8)+floor` recebe
`noncontracting_refinement`; produto dos incrementos menor que `-floor²`
recebe `oscillatory_refinement`. Esses casos não autorizam uma classificação
resolvida. O `1e-8` evita classificar arredondamento relativo como crescimento;
não modifica gates do solver ou tolerâncias geométricas.

Para o par A/B, as diferenças `B-A` são conservadas nos **dois níveis mais
finos, em cada método**. Os padrões são:

| Classificação | Critério |
|---|---|
| `numerically_distinguishable` | Histórias completas, sem limitações ou rótulos de refinamento; todas as diferenças têm o mesmo sinal; o menor módulo excede `max(factor*(U(A)+U(B)), floor)`. |
| `inconsistent_across_refinement_or_methods` | Há diferenças maiores que `floor` com sinais opostos entre métodos ou entre os dois níveis mais finos. |
| `numerically_unresolved` | Evidência ausente, inválida, incompatível ou limitada; ou separação insuficiente segundo o critério acima. |

Dados incompletos/limitados têm precedência e impedem inferir até uma troca de
ordenação. Em histórias completas, uma troca observada recebe a classificação
inconsistente mesmo quando também existe refinamento oscilatório. O padrão
`RankingPolicy(safety_factor=2, comparison_floor=1e-12)` foi declarado antes
da campanha. Esses parâmetros são escolhas heurísticas explícitas e podem ser
passados às duas funções; não foram calibrados para aprovar candidatos.

O resultado inclui `relation` (`a_greater`/`b_greater` apenas se resolvido),
diferenças, indicadores, margem, política e rótulos. `to_dict()` permite
serializar as histórias. Identificadores são rótulos escolhidos pelo chamador;
não substituem o ID científico do runner.

## Alcance e limites

Concordância entre métodos pode conter erro comum. Refinamento finito pode
preceder oscilações ou outro regime assintótico. O indicador não limita erro
de geometria, modelo, pressão local, Reynolds, viscosidade, separação ou estol.
Mesmo uma ordenação estável e numericamente distinguível não estabelece
desempenho físico nem autoriza otimização automática. Um benchmark isolado
não transforma U em estimativa universal. O campo
`indicator_kind="heuristic_not_error_bound"` conserva essa distinção nos dados.

O harness `examples/verify_ranking_reliability.py` realiza estudos somente
quando solicitado. A campanha e suas limitações estão em `VALIDATION_V032.md`.
O `PanelSolver` continua executável sem importar este módulo e sem refinamento
automático ou custo adicional por execução.
