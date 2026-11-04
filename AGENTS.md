# Instruções permanentes — AeroOpt

- Consulte `docs/PROJECT_CONTEXT.md` antes de alterar a arquitetura.
- Preserve e consulte `docs/GEOMETRY_CONVENTIONS.md` para alterações geométricas.
- Use exclusivamente `.venv\Scripts\python.exe` para operações Python do projeto.
  Não recrie a `.venv` nem instale dependências globalmente.
- Execute os testes após modificações: `.\.venv\Scripts\python.exe -m pytest -q`.
- Atualize a documentação quando decisões técnicas ou a API forem alteradas.
- Não confunda resultados geométricos com resultados aerodinâmicos.
- Não modifique dados de entrada silenciosamente. Não suavize, interpole, reordene
  ou redistribua pontos automaticamente.
- Não introduza funcionalidades fora do escopo sem justificativa documentada.
- Não execute comandos destrutivos sem autorização. Preserve o histórico Git,
  não altere remotos e não crie commits automaticamente.
- Mantenha arquivos gerados, caches e temporários dentro do workspace.
- A v0.1 é exclusivamente geométrica. CFD, solvers aerodinâmicos, NSGA-II e
  modelos substitutos pertencem a versões futuras.

