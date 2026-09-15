# 发布与出版（GitHub → Zenodo DOI）

本目录是发布用**归档树**：只含代码与说明文档，不含数据与图件。

## 一、更新 GitHub 仓库

```bash
# 在归档目录内初始化 / 关联远端（首次）
cd zenodo_archive/dlbcl-hmgb1-havcr2-mafb-analysis-v1.0.0
git init -b main
git remote add origin https://github.com/tjwu87/dlbcl-hmgb1-havcr2-mafb-analysis.git

git add -A
git commit -m "v1.0.0: reorganise by paper section (Results 1-6); refresh analysis code"
git push -u origin main --force      # 覆盖远端默认分支，抹掉旧版内容
```

> 若希望旧代码在 GitHub 上**完全不可见**（含历史提交），最干净的做法是
> 删除仓库后以同名重建，再 push 本目录。

## 二、保持同一个 DOI 更新 Zenodo（不新增版本）

Zenodo 的文件只能在**发布后 30 天内**自助替换：

1. 打开记录页 → **Edit** → 展开 **Edit files** → **Edit published files**
2. 替换文件（上传新的 ZIP，删除旧 ZIP）
3. 点击 **Publish**

官方说明：*"Publishing the draft will not change the DOI."*
（草稿须在原始发布后 45 天内发布。）

⚠️ **不要再创建新的 GitHub Release**：本记录由 GitHub Release 集成建立，
每发一个新 Release，Zenodo 都会自动生成**新版本 + 新 DOI**。
更新仓库默认分支不会影响 Zenodo。

若 "Edit published files" 不可用（GitHub 关联记录可能受限），退路是：
在记录页把旧文件设为**限制访问**（restrict public access），
使 DOI 记录仍在、元数据公开，但旧代码无法下载；必要时联系 Zenodo support。

## 三、若将来需要正式发布新版本

Zenodo 记录页 → **New version** → 上传本 ZIP → Publish。
会铸造新的 **Version DOI**；原 **Concept DOI** 不变、并始终指向最新版。

## 四、归档完整性

见 `MANIFEST.md`（逐文件 SHA256）。
